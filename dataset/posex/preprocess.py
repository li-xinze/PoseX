import os
import json
import inspect
import requests
import subprocess
import pandas as pd
from pymol import cmd
from typing import Any
from rdkit import Chem
from functools import partial
from multiprocessing import Pool
from collections import defaultdict
from posex.utils import bcif2cif, my_tqdm
from pdbecif.mmcif_io import CifFileReader
from posex.utils import run_in_tmp_dir, get_ccd_instance_map


NUM_CPUS = 100
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class DataPreprocessor():
    """
    A class for downloading and preprocessing data

    Args:
        pdbid_list : list of pdbids
    """
    def __init__(self, pdbid_list: list[str], **kwargs: Any):
        self.bcif_dir = kwargs.get("bcif_dir")
        self.cif_dir = kwargs.get("cif_dir")
        self.vs_dir = kwargs.get("vs_dir")
        self.ccd_dir = kwargs.get("ccd_dir")
        self.ccd_path = kwargs.get("ccd_path")
        self.lig_dir = kwargs.get("lig_dir")
        self.molecule_dir = kwargs.get("molecule_dir")
        self.components_path = kwargs.get("components_path")
        self.pdbid_list = pdbid_list
        
    def _get_func_name(self, idx=1) -> str:
        """return the name of a funciton in the current call stack
        """
        return inspect.stack()[idx].function

    def download_bcif(self) -> None:
        os.makedirs(self.bcif_dir, exist_ok=True)
        for pdbid in my_tqdm(self.pdbid_list, desc=self._get_func_name()):
            if os.path.exists(f"{self.bcif_dir}/{pdbid}.bcif"):
                continue
            try:
                url = f"https://models.rcsb.org/{pdbid}.bcif"
                subprocess.run(["wget", url], cwd=self.bcif_dir, check=True)
            except Exception as e:
                print(pdbid, e)

    def convert_bif_to_cif(self, num_cpus=NUM_CPUS) -> None:
        os.makedirs(self.cif_dir, exist_ok=True)
        bcif2cif_ = partial(bcif2cif, bcif_dir=self.bcif_dir, cif_dir=self.cif_dir)
        with Pool(processes=num_cpus) as pool:
            for bcif_file, success in my_tqdm(pool.imap_unordered(bcif2cif_, os.listdir(self.bcif_dir)), desc=self._get_func_name()):
                if not success:
                    print(f" {bcif_file} failed")
              
    def download_validation_score(self) -> None:
        os.makedirs(self.vs_dir, exist_ok=True)
        url = 'https://data.rcsb.org/graphql'
        with open(f"{PROJECT_DIR}/template/vs_query.txt", "r") as f:
            lines = f.readlines()
            query_string = "".join([line.strip() for line in lines])
        for pdbid in my_tqdm(self.pdbid_list, desc=self._get_func_name()):
            json_path = f"{self.vs_dir}/{pdbid}.json"
            if os.path.exists(json_path):
                continue
            content = {
                "query": query_string,
                "variables": {"id": pdbid}
            }
            response = requests.post(url, json=content)
            with open(json_path, 'w') as f:
                json.dump(response.json(), f)

    def query_ccd(self, ccd: str) -> tuple[float, str]:
        """Query formula_weight and inchi of a ligand using rcsb API

        Args:
            ccd (str): chemical descriptions of a ligand

        Returns:
            tuple[float, str]: formula_weight and inchi of the ligand
        """
        json_path = os.path.join(self.ccd_dir, f"{ccd}.json")
        if not os.path.exists(json_path):
            url = 'https://data.rcsb.org/graphql'
            with open(f"{PROJECT_DIR}/template/ccd_query.txt", "r") as f:
                lines = f.readlines()
            query_string = "".join([line.strip() for line in lines])
            content = {
                "query": query_string,
                "variables": {"id": ccd}
            }
            response = requests.post(url, json=content)
            with open(json_path, 'w') as f:
                json.dump(response.json(), f)

        with open(json_path, "r") as f:
            info = json.load(f)
            formula_weight = info["data"]["chem_comp"]["chem_comp"]["formula_weight"]
            inchi = None
            pdbx_chem_comp_descriptor = info["data"]["chem_comp"]["pdbx_chem_comp_descriptor"]
            if pdbx_chem_comp_descriptor is not None:
                for item in pdbx_chem_comp_descriptor:
                    if item["type"] == "InChI":
                        inchi = item["descriptor"]
                        break
        return formula_weight, inchi

    def _download_ccd_dict(self) -> str:
        """Download components.cif
        """
        os.makedirs(self.ccd_dir, exist_ok=True)
        if not os.path.exists(self.components_path):
            print("downloading components.cif...")
            url = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif"
            subprocess.run(["wget", url], cwd=self.ccd_dir, check=True)
            print("finish")
        print("finish")

    def build_ccd_table(self) -> None:
        """Create a ccd table and save it to self.ccd_path
        """
        if not os.path.exists(self.ccd_path):
            self._download_ccd_dict()
            cfr = CifFileReader()
            cif_obj = cfr.read(self.components_path, ignore = ['_atom_site'])
            df_dict = defaultdict(list)
            for ccd, values in cif_obj.items():
                if ccd == "UNL": continue
                descriptor_types = values["_pdbx_chem_comp_descriptor"]["type"]
                inchi_idx = -1
                for idx, dtype in enumerate(descriptor_types):
                    if dtype == "InChI":
                        inchi_idx = idx
                        break
                if inchi_idx == -1:
                    print(f"no InChI for {ccd}")
                else:
                    inchi = values["_pdbx_chem_comp_descriptor"]["descriptor"][inchi_idx]
                    df_dict["CCD"].append(ccd)
                    df_dict["InChI"].append(inchi)
                    df_dict["MOLWT"].append(values["_chem_comp"]["formula_weight"])
        else:
            df_dict = pd.read_csv(self.ccd_path).to_dict(orient="list")
            df_dict["CCD"] = [ccd.strip("''") for ccd in df_dict["CCD"]]
        
        # get all ccds
        ccd_set = set()
        for pdbid in self.pdbid_list:
            json_path = os.path.join(self.vs_dir, f"{pdbid}.json")
            with open(json_path, "r") as f:
                info = json.load(f)
                try:
                    for nonpolymer_entity in info["data"]["entry"]["nonpolymer_entities"]:
                        ccd = nonpolymer_entity["nonpolymer_comp"]["chem_comp"]["id"]
                        ccd_set.add(ccd)
                except:
                    print(json_path)
                    raise

        uncover_ccds = ccd_set - set(df_dict["CCD"])
        for ccd in my_tqdm(uncover_ccds, desc=self._get_func_name()):
            formula_weight, inchi = self.query_ccd(ccd)
            if formula_weight is not None and inchi is not None:
                df_dict["CCD"].append(ccd)
                df_dict["InChI"].append(inchi)
                df_dict["MOLWT"].append(formula_weight)
        df = pd.DataFrame(df_dict)
        
        df.CCD = df.CCD.apply(lambda x: f"\'{x}\'")
        df.to_csv(self.ccd_path, na_rep=None, index=False)

    @run_in_tmp_dir
    def download_ligand(self, num_cpus=20) -> dict[str, dict[str, list[str]]]:
        """Download sdf files from rcsb

        Returns:
            dict[str, dict[str, list[str]]]:
                dict of pdb to dict of ccd to list of asym_ids
        """
        pdb_ccd_instance_map = {}
        fetch_ligands_ = partial(fetch_ligands, cif_dir=self.cif_dir, lig_dir=self.lig_dir)
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(fetch_ligands_, self.pdbid_list)
            for pdb, ccd_instance_map in my_tqdm(pool_iter, total=len(self.pdbid_list), desc="downloading ligand"):
                pdb_ccd_instance_map[pdb] = ccd_instance_map
        return pdb_ccd_instance_map
        
    @run_in_tmp_dir
    def extract_molecule(self) -> None:
        """Extract sdf files of organic_molecule and metal_ion from cif files
        """
        os.makedirs(self.molecule_dir, exist_ok=True)
        for pdb in my_tqdm(self.pdbid_list, desc=self._get_func_name()):
            cif_path = os.path.join(self.cif_dir, f'{pdb}.cif')
            organic_molecule_path = f"{self.molecule_dir}/{pdb}_organic_molecule.cif"
            metal_ion_path = f"{self.molecule_dir}/{pdb}_metal_ion.cif"
            if os.path.exists(organic_molecule_path) and os.path.exists(metal_ion_path):
                continue
            cmd.load(cif_path, pdb)
            cmd.select("organic_molecule", "organic")
            cmd.select("metal_ion", "metal")
            cmd.save(organic_molecule_path, "organic_molecule")
            cmd.save(metal_ion_path, "metal_ion")
            cmd.delete("all")

    def run(self) -> dict[str, dict[str, list[str]]]:
        """Preprocess

        Returns:
            dict[str, dict[str, list[str]]]:
                dict of pdb to dict of ccd to list of asym_ids
        """
        self.download_bcif()
        self.convert_bif_to_cif()
        self.download_validation_score()
        self.build_ccd_table()
        self.extract_molecule()
        pdb_ccd_instance_map = self.download_ligand()
        return pdb_ccd_instance_map

def fetch_ligands(pdb, cif_dir, lig_dir):
    _, ccd_instance_map = get_ccd_instance_map(pdb, cif_dir)
    for ccd, asym_ids in ccd_instance_map.items():
        mols = []
        ligand_path = f"{lig_dir}/{pdb}_{ccd}.sdf"
        if os.path.exists(ligand_path):
            continue
        for asym_id in asym_ids:
            instance_name = f"{pdb}_{ccd}_{asym_id}"
            tmp_ligand_path = f"{lig_dir}/{instance_name}.sdf"
            url = f"https://models.rcsb.org/v1/{pdb}/ligand?label_comp_id={ccd}&label_asym_id={asym_id}&encoding=sdf"
            try:
                r = requests.get(url)
                open(tmp_ligand_path , 'wb').write(r.content)
            except Exception as e:
                print(e, pdb, ccd, asym_id)
            mol = Chem.SDMolSupplier(tmp_ligand_path)[0]
            os.remove(tmp_ligand_path)
            if mol is not None:
                params = Chem.RemoveHsParameters()
                params.removeDegreeZero = True
                mol = Chem.RemoveHs(mol, params)
                props = list(mol.GetPropNames())
                for prop in props:
                    mol.ClearProp(prop)
                mol.SetProp('_Name', instance_name)
                mols.append(mol)
            else:
                # All ligand SDF files can be loaded with RDKit and pass its sanitization
                mols = []
                break
        w = Chem.SDWriter(ligand_path)
        for mol in mols:
            w.write(mol)
        w.close()
    return pdb, ccd_instance_map