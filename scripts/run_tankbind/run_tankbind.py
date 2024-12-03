import os
import sys
import torch
import shutil
import logging
import argparse
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from Bio.PDB import PDBParser
from torch_geometric.loader import DataLoader

torch.set_num_threads(1)

class TankBindRunner():
    def __init__(self, args) -> None:
        self.input_dir = args.input_dir
        self.output_dir = args.output_dir
        self.tankbind_exec_dir = args.tankbind_exec_dir
        self.p2rank_exec_path = args.p2rank_exec_path
        input_data_path = os.path.join(self.input_dir, "data.csv")
        self.input_data = pd.read_csv(input_data_path)
        self.pdb_list = self.input_data.PDB_CCD_ID.values
        self.rdkit_folder = f"{self.output_dir}/rdkit"
        os.makedirs(self.rdkit_folder, exist_ok=True)

    def _create_data(self, pockets_dict, protein_dict):
        info = []
        for pdb in self.pdb_list:
            protein_name = pdb
            compound_name = pdb
            pocket = pockets_dict[pdb].head(10)
            pocket.columns = pocket.columns.str.strip()
            pocket_coms = pocket[['center_x', 'center_y', 'center_z']].values
            # native block.
            info.append([protein_name, compound_name, pdb, None, True, False])
            # protein center as a block.
            protein_com = protein_dict[protein_name][0].numpy().mean(axis=0).astype(float).reshape(1, 3)
            info.append([protein_name, compound_name, pdb+"_c", protein_com, 0, False, False])
            for idx, pocket_line in pocket.iterrows():
                pdb_idx = f"{pdb}_{idx}"
                info.append([protein_name, compound_name, pdb_idx, pocket_coms[idx].reshape(1, 3), False, False])
        info = pd.DataFrame(info, columns=['protein_name', 'compound_name', 'pdb', 'pocket_com', 'affinity',
                                        'use_compound_com', 'use_whole_protein'])
        return info

    def _predict_protein_feature(self):
        protein_dict = {}
        for pdb in self.pdb_list:
            proteinFile = f"{self.input_dir}/protein_remove_extra_chains_10A/{pdb}_protein.pdb"
            parser = PDBParser(QUIET=True)
            s = parser.get_structure(pdb, proteinFile)
            res_list = get_clean_res_list(s.get_residues(), verbose=False, ensure_ca_exist=True)
            protein_dict[pdb] = get_protein_feature(res_list)
        return protein_dict

    def _predict_ligand_feature(self):
        compound_dict = {}
        for _, row in self.input_data.iterrows():
            pdb = row["PDB_CCD_ID"]
            sdf_path = row["LIGAND_SDF_PATH"]
            mol = Chem.SDMolSupplier(sdf_path)[0]
            smiles = Chem.MolToSmiles(mol)
            rdkit_mol_path = f"{self.rdkit_folder}/{pdb}_ligand.sdf"
            generate_sdf_from_smiles_using_rdkit(smiles, rdkit_mol_path, shift_dis=0)
            mol = Chem.SDMolSupplier(rdkit_mol_path)[0]
            compound_dict[pdb] = extract_torchdrug_feature_from_mol(mol, has_LAS_mask=True)  # self-dock set has_LAS_mask to true
        return compound_dict

    def _predict_pockets(self):
        # predict pockets by p2rank 
        p2rank_prediction_folder = f"{self.input_dir}/p2rank_protein_remove_extra_chains_10A"
        os.system(f"mkdir -p {p2rank_prediction_folder}")
        ds = f"{p2rank_prediction_folder}/protein_list.ds"
        with open(ds, "w") as out:
            for pdb in self.pdb_list:
                out.write(f"../protein_remove_extra_chains_10A/{pdb}_protein.pdb\n")
        cmd = ["bash", self.p2rank_exec_path, "predict", ds, "-o", f"{p2rank_prediction_folder}/p2rank", "-threads", "8"]
        subprocess.run(cmd, check=True)

        # handle predictions
        d_list = []
        for name in self.pdb_list:
            p2rankFile = f"{self.input_dir}/p2rank_protein_remove_extra_chains_10A/p2rank/{name}_protein.pdb_predictions.csv"
            d = pd.read_csv(p2rankFile)
            d.columns = d.columns.str.strip()
            d_list.append(d.assign(name=name))
        d = pd.concat(d_list).reset_index(drop=True)
        d.reset_index(drop=True).to_feather(f"{self.input_dir}/p2rank_result.feather")
        d = pd.read_feather(f"{self.input_dir}/p2rank_result.feather")

        pockets_dict = {}
        for name in self.pdb_list:
            pockets_dict[name] = d[d.name == name].reset_index(drop=True)
        return pockets_dict

    def _save_conformation(self, dataset, chosen, y_pred_list):
        device = "cpu"
        for _, line in tqdm(chosen.iterrows(), total=chosen.shape[0]):
            name = line['compound_name']
            dataset_index = line['dataset_index']
            coords = dataset[dataset_index].coords.to(device)
            protein_nodes_xyz = dataset[dataset_index].node_xyz.to(device)
            n_compound = coords.shape[0]
            n_protein = protein_nodes_xyz.shape[0]
            y_pred = y_pred_list[dataset_index].reshape(n_protein, n_compound).to(device)
            y = dataset[dataset_index].dis_map.reshape(n_protein, n_compound).to(device)
            compound_pair_dis_constraint = torch.cdist(coords, coords)
            rdkit_mol_path = f"{self.rdkit_folder}/{name}_ligand.sdf"
            mol = Chem.SDMolSupplier(rdkit_mol_path)[0]
            LAS_distance_constraint_mask = get_LAS_distance_constraint_mask(mol).bool()
            pred_dist_info = get_info_pred_distance(coords, y_pred, protein_nodes_xyz, compound_pair_dis_constraint,
                                        LAS_distance_constraint_mask=LAS_distance_constraint_mask,
                                        n_repeat=1, show_progress=False)
            toFile = f'{self.output_dir}/{name}_tankbind_chosen.sdf'
            new_coords = pred_dist_info.sort_values("loss")['coords'].iloc[0].astype(np.double)
            write_with_new_coords(mol, new_coords, toFile)

    def process_testset(self):
        toFolder = f"{self.input_dir}/protein_remove_extra_chains_10A/"
        os.makedirs(toFolder, exist_ok=True)
        for _, row in self.input_data.iterrows():
            cutoff = 10
            itemname = row["PDB_CCD_ID"]
            toFile = f"{toFolder}/{itemname}_protein.pdb"
            x = (row["PROTEIN_PDB_PATH"], row["LIGAND_SDF_PATH"], cutoff, toFile)
            select_chain_within_cutoff_to_ligand_v2(x)
        pockets_dict = self._predict_pockets()
        protein_dict = self._predict_protein_feature()
        compound_dict = self._predict_ligand_feature()
        data = self._create_data(pockets_dict, protein_dict)
        dataset_dir = f"{self.input_dir}/dataset"
        if os.path.exists(dataset_dir):
            shutil.rmtree(dataset_dir)
        os.makedirs(dataset_dir)
        testset = TankBindDataSet(dataset_dir, data=data, protein_dict=protein_dict, compound_dict=compound_dict)
        testset = TankBindDataSet(dataset_dir, proteinMode=0, compoundMode=1, pocket_radius=20, predDis=True)
        return testset

    def predict(self, dataset, device="cpu"):
        data_loader = DataLoader(dataset, batch_size=1, 
                                follow_batch=['x', 'y', 'compound_pair'], shuffle=False, num_workers=8, pin_memory=True)
        logging.basicConfig(level=logging.INFO)
        model = get_model(0, logging, device)
        model.eval()
        model.load_state_dict(torch.load(f"{self.tankbind_exec_dir}/saved_models/self_dock.pt", map_location=device))
        y_pred_list, affinity_pred_list = [], []
        for data in tqdm(data_loader):
            data = data.to(device)
            with torch.no_grad():
                y_pred, affinity_pred = model(data)
            affinity_pred_list.append(affinity_pred.detach().cpu())
            for i in range(data.y_batch.max() + 1):
                y_pred_list.append((y_pred[data['y_batch'] == i]).detach().cpu())
        affinity_pred_list = torch.cat(affinity_pred_list)
        output_info_chosen = dataset.data
        output_info_chosen['affinity'] = affinity_pred_list
        output_info_chosen['dataset_index'] = range(len(output_info_chosen))
        output_info_chosen = output_info_chosen.query("not use_compound_com").reset_index(drop=True)
        chosen = output_info_chosen.loc[output_info_chosen.groupby(['protein_name', 'compound_name'], sort=False)['affinity'].agg('idxmax')].reset_index()
        self._save_conformation(dataset, chosen, y_pred_list)

    def run(self):
        testset = self.process_testset()
        self.predict(testset)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tankbind_exec_dir", type=str, required=True, help="Path to the TankBind codebase")
    parser.add_argument("--p2rank_exec_path", type=str, required=True, help="Path to the p2rank_exec_path")
    parser.add_argument("--input_dir", type=str, required=True, help="Path to the input dir")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--gpu_id", type=int, required=True, help="GPU ID")
    args = parser.parse_args()

    tankbind_module_dir = os.path.join(args.tankbind_exec_dir, "tankbind")
    sys.path.insert(0, tankbind_module_dir)
    from data import TankBindDataSet
    from model import get_model
    from feature_utils import select_chain_within_cutoff_to_ligand_v2, get_protein_feature, \
        get_clean_res_list, extract_torchdrug_feature_from_mol, generate_sdf_from_smiles_using_rdkit
    from generation_utils import get_LAS_distance_constraint_mask, get_info_pred_distance, write_with_new_coords

    tankbind_runner = TankBindRunner(args)
    tankbind_runner.run()
