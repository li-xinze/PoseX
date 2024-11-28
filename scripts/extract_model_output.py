import argparse
import os

import numpy as np
import pandas as pd
import prody
from Bio.PDB import MMCIFParser, PDBIO
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem


def extract_alphafold3_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Posebusters Data: ", len(docking_data))

    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"].lower()
        ligand_smiles = row["LIGAND_SMILES"]

        if not os.path.exists(os.path.join(args.output_folder, f"{pdb_ccd_id}")):
            print(f"Directory {pdb_ccd_id} does not exist")
            continue

        # Parse the MMCIFFile
        pdb = prody.parseMMCIF(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.cif"))
        protein = pdb.select("protein")
        ligand = pdb.select("not (protein or nucleotide or water)")

        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_protein.pdb"), protein)
        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), ligand)

        # Replace the ligand ID from LIG_Z to LIG
        with open(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), 'r') as f:
            content = f.read()
        content = content.replace("LIG_Z", "LIG")
        with open(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), 'w') as f:
            f.write(content)

        # Convert the ligand PDB to SDF
        mol = Chem.MolFromPDBFile(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), sanitize=False)
        template = AllChem.MolFromSmiles(ligand_smiles)
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        writer = Chem.SDWriter(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.sdf"))
        writer.write(mol)
        writer.close()


def extract_chai_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Posebusters Data: ", len(docking_data))

    error_process_ligand_ids = []

    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        ligand_smiles = row["LIGAND_SMILES"]

        if not os.path.exists(os.path.join(args.output_folder, f"{pdb_ccd_id}")):
            print(f"Directory {pdb_ccd_id} does not exist")
            continue

        # Find the best model
        model_scores = []
        for model_idx in range(5):
            score_file = os.path.join(args.output_folder, f"{pdb_ccd_id}/scores.model_idx_{model_idx}.npz")
            model_scores.append(np.load(score_file)["aggregate_score"])
        best_model_idx = np.argmax(model_scores)

        # Parse the MMCIFFile to PDB
        parser = MMCIFParser()
        structure = parser.get_structure(pdb_ccd_id, os.path.join(args.output_folder, f"{pdb_ccd_id}/pred.model_idx_{best_model_idx}.cif"))
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))

        # Convert the ligand atoms to HETATM
        pdb = PandasPdb().read_pdb(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))
        ligand_atoms = pdb.df["ATOM"][pdb.df["ATOM"]["atom_name"] == "LIG"]
        ligand_indices = ligand_atoms.index
        pdb.df["ATOM"] = pdb.df["ATOM"].drop(ligand_indices)
        ligand_atoms.record_name = "HETATM"
        pdb.df["HETATM"] = ligand_atoms
        pdb.to_pdb(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))
        
        # Parse the PDBFile
        pdb = prody.parsePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))
        protein = pdb.select("protein")
        ligand = pdb.select("not (protein or nucleotide or water)")

        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_protein.pdb"), protein)
        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), ligand)

        # Convert the ligand PDB to SDF
        mol = Chem.MolFromPDBFile(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), sanitize=False)
        try:
            template = AllChem.MolFromSmiles(ligand_smiles)
            mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        except Exception as e:
            mol = Chem.MolFromPDBFile(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), sanitize=False)
            error_process_ligand_ids.append(pdb_ccd_id)
            print(f"Error processing ligand for {pdb_ccd_id}: {e}")
        writer = Chem.SDWriter(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.sdf"))
        writer.write(mol)
        writer.close()

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def main(args: argparse.Namespace):
    if args.model_type == "alphafold3":
        extract_alphafold3_output(args)
    elif args.model_type == "chai":
        extract_chai_output(args)
    else:
        raise ValueError(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the alphafold3 output folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()
    
    main(args)
