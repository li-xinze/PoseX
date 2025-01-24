import argparse
import os
import glob
import shutil
import numpy as np
import pandas as pd
import prody
from Bio.PDB import MMCIFParser, PDBIO
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem



def convert_ligand_pdb_to_sdf(model_output_folder: str, pdb_ccd_id: str, ligand_smiles: str) -> bool:
    """Convert the ligand PDB to SDF

    Args:
        model_output_folder (str): Path to the model output folder
        pdb_ccd_id (str): PDB_CCD_ID
        ligand_smiles (str): Ligand SMILES

    Returns:
        bool: Whether the conversion is successful
    """
    mol = Chem.MolFromPDBFile(os.path.join(model_output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), sanitize=False)
    try:
        template = AllChem.MolFromSmiles(ligand_smiles)
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    except Exception as e:
        mol = Chem.MolFromPDBFile(os.path.join(model_output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), sanitize=False)
        print(f"Error assigning bond orders for ligand {pdb_ccd_id}: {e}")
    
    try:
        test_mol = Chem.Mol(mol)
        Chem.SanitizeMol(test_mol)
    except Exception as e:
        print(f"Sanitization failed for ligand {pdb_ccd_id}: {e}")
        return False

    writer = Chem.SDWriter(os.path.join(model_output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.sdf"))
    writer.write(mol)
    writer.close()
    return True


def extract_alphafold3_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))

    error_process_ligand_ids = []

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
        convert_success = convert_ligand_pdb_to_sdf(args.output_folder, pdb_ccd_id, ligand_smiles)
        if not convert_success:
            print(f"Error processing ligand for {pdb_ccd_id}")
            error_process_ligand_ids.append(pdb_ccd_id)
            continue

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_chai_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))

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
        print("best_model_idx: ", best_model_idx)

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname in [f"LIG{i}" for i in range(1, 10)]:
                        residue.resname = "LIG"
                    for atom in residue:
                        if "_" in atom.name:
                            atom.name = atom.name.split("_")[0]
                            atom.fullname = atom.name

        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))
        
        # Parse the PDBFile
        pdb = prody.parsePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))
        protein = pdb.select("protein")
        ligand = pdb.select("not (protein or nucleotide or water)")

        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_protein.pdb"), protein)
        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), ligand)

        # Convert the ligand PDB to SDF
        convert_success = convert_ligand_pdb_to_sdf(args.output_folder, pdb_ccd_id, ligand_smiles)
        if not convert_success:
            print(f"Error processing ligand for {pdb_ccd_id}")
            error_process_ligand_ids.append(pdb_ccd_id)
            continue

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_boltz_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))

    error_process_ligand_ids = []

    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        ligand_smiles = row["LIGAND_SMILES"]

        if not os.path.exists(os.path.join(args.output_folder, f"{pdb_ccd_id}")):
            print(f"Directory {pdb_ccd_id} does not exist")
            continue

        cif_path = os.path.join(args.output_folder, f"{pdb_ccd_id}/boltz_results_{pdb_ccd_id}/predictions/{pdb_ccd_id}/{pdb_ccd_id}_model_0.cif")
        if not os.path.exists(cif_path):
            print(f"CIF file {cif_path} does not exist")
            continue

        parser = MMCIFParser()
        structure = parser.get_structure(pdb_ccd_id, cif_path)
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))

        # Parse the PDBFile
        pdb = prody.parsePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model.pdb"))
        protein = pdb.select("protein")
        ligand = pdb.select("not (protein or nucleotide or water)")

        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_protein.pdb"), protein)
        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), ligand)

        # Convert the ligand PDB to SDF
        convert_success = convert_ligand_pdb_to_sdf(args.output_folder, pdb_ccd_id, ligand_smiles)
        if not convert_success:
            print(f"Error processing ligand for {pdb_ccd_id}")
            error_process_ligand_ids.append(pdb_ccd_id)
            continue

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_rfaa_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))

    error_process_ligand_ids = []
    
    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        ligand_smiles = row["LIGAND_SMILES"]

        if not os.path.exists(os.path.join(args.output_folder, f"{pdb_ccd_id}")):
            print(f"Directory {pdb_ccd_id} does not exist")
            continue

        pdb_path = os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}.pdb")
        if not os.path.exists(pdb_path):
            print(f"PDB file {pdb_path} does not exist")
            continue

        # Parse the PDBFile
        pdb = prody.parsePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}.pdb"))
        protein = pdb.select("protein")
        ligand = pdb.select("not (protein or nucleotide or water)")

        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_protein.pdb"), protein)
        prody.writePDB(os.path.join(args.output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand.pdb"), ligand)

        # Convert the ligand PDB to SDF
        convert_success = convert_ligand_pdb_to_sdf(args.output_folder, pdb_ccd_id, ligand_smiles)
        if not convert_success:
            print(f"Error processing ligand for {pdb_ccd_id}")
            error_process_ligand_ids.append(pdb_ccd_id)
            continue

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_tankbind_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))
    error_process_ligand_ids = []
    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        output_dir = os.path.join(args.output_folder, pdb_ccd_id)
        os.makedirs(output_dir, exist_ok=True)
        output_pdb_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_protein.pdb")
        output_sdf_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_ligand.sdf")
        try:
            shutil.copy(row["PROTEIN_PDB_PATH"], output_pdb_path)
            shutil.copy(f"{args.output_folder}/{pdb_ccd_id}_tankbind_chosen.sdf", output_sdf_path)
        except Exception as e:
            error_process_ligand_ids.append(pdb_ccd_id)
            print(f"Error processing ligand for {pdb_ccd_id}: {e}")   

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_dynamicbind_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))
    error_process_ligand_ids = []
    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        output_dir = os.path.join(args.output_folder, pdb_ccd_id)
        output_pdb_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_protein.pdb")
        output_sdf_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_ligand.sdf")
        try:
            matching_protein_files = glob.glob(os.path.join(args.output_folder, pdb_ccd_id, "index0_idx_0", f"rank1_receptor*"))
            matching_ligand_files = glob.glob(os.path.join(args.output_folder, pdb_ccd_id, "index0_idx_0", f"rank1_ligand*"))
            assert len(matching_protein_files) == 1
            assert len(matching_ligand_files) == 1
            shutil.copy(matching_protein_files[0], output_pdb_path)
            shutil.copy(matching_ligand_files[0], output_sdf_path)
        except Exception as e:
            error_process_ligand_ids.append(pdb_ccd_id)
            print(f"Error processing ligand for {pdb_ccd_id}: {e}")   

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_diffdock_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))
    error_process_ligand_ids = []
    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        output_dir = os.path.join(args.output_folder, pdb_ccd_id)
        os.makedirs(output_dir, exist_ok=True)
        output_pdb_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_protein.pdb")
        output_sdf_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_ligand.sdf")
        try:
            shutil.copy(row["PROTEIN_PDB_PATH"], output_pdb_path)
            shutil.copy(f"{args.output_folder}/{pdb_ccd_id}/rank1.sdf", output_sdf_path)
        except Exception as e:
            error_process_ligand_ids.append(pdb_ccd_id)
            print(f"Error processing ligand for {pdb_ccd_id}: {e}")   

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def extract_fabind_output(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    print("Number of Benchmark Data: ", len(docking_data))
    error_process_ligand_ids = []
    for _, row in docking_data.iterrows():
        print(f"Processing {row['PDB_CCD_ID']}")
        pdb_ccd_id = row["PDB_CCD_ID"]
        output_dir = os.path.join(args.output_folder, pdb_ccd_id)
        os.makedirs(output_dir, exist_ok=True)
        output_pdb_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_protein.pdb")
        output_sdf_path = os.path.join(output_dir, f"{pdb_ccd_id}_model_ligand.sdf")
        shutil.copy(row["PROTEIN_PDB_PATH"], output_pdb_path)
        try:
            matching_files = glob.glob(os.path.join(args.output_folder, f'{pdb_ccd_id}_*'))
            assert len(matching_files) == 1
            shutil.copy(matching_files[0], output_sdf_path)
        except Exception as e:
            error_process_ligand_ids.append(pdb_ccd_id)
            print(f"Error processing ligand for {pdb_ccd_id}: {e}")   

    print(f"Number of Error Process Ligand IDs: {len(error_process_ligand_ids)}")
    print(f"Error Process Ligand IDs: {error_process_ligand_ids}")


def main(args: argparse.Namespace):
    if args.model_type == "alphafold3":
        extract_alphafold3_output(args)
    elif args.model_type == "chai":
        extract_chai_output(args)
    elif args.model_type == "boltz":
        extract_boltz_output(args)
    elif args.model_type == "rfaa":
        extract_rfaa_output(args)
    elif args.model_type == "tankbind":
        extract_tankbind_output(args)
    elif args.model_type in ["diffdock", "diffdock_l"]:
        extract_diffdock_output(args)
    elif args.model_type == "fabind":
        extract_fabind_output(args)
    elif args.model_type == "dynamicbind":
        extract_dynamicbind_output(args)
    else:
        raise ValueError(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the model output folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()
    
    main(args)
