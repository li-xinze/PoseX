import argparse
import json
import os
import shutil
import string

import pandas as pd
import yaml
from tqdm.auto import tqdm


def generate_alphafold3_input(args: argparse.Namespace):
    """Generate AlphaFold3 input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        task_dict = {"name": task_name, "sequences": []}
        
        protein_sequences = row["PROTEIN_SEQUENCE"].split("|")
        # Skip if there are more than 20 chains
        if len(protein_sequences) > 20:
            print(f"Skipping {row['PDB_CCD_ID']} because it has more than 20 chains")
            continue
        
        # Add the protein sequences
        for sequence_idx, protein_sequence in enumerate(protein_sequences):
            sequence_dict = {
                "protein": {
                    "id": string.ascii_uppercase[sequence_idx],
                    "sequence": protein_sequence
                }
            }
            task_dict["sequences"].append(sequence_dict)

        # Add the ligand
        task_dict["sequences"].append({
            "ligand": {
                "id": "Z",
                "smiles": row["LIGAND_SMILES"]
            }
        })

        task_dict["modelSeeds"] = [42]
        task_dict["dialect"] = "alphafold3"
        task_dict["version"] = 1

        part_folder = os.path.join(args.output_folder, f"part_{data_idx % 4}")
        os.makedirs(part_folder, exist_ok=True)
        with open(os.path.join(part_folder, f"{task_name}.json"), "w") as f:
            json.dump(task_dict, f, indent=2)


def generate_chai_input(args: argparse.Namespace):
    """Generate CHAI input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    for _, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        fasta_string = ""

        protein_sequences = row["PROTEIN_SEQUENCE"].split("|")
        # Skip if there are more than 20 chains
        if len(protein_sequences) > 20:
            print(f"Skipping {row['PDB_CCD_ID']} because it has more than 20 chains")
            continue

        # Add the protein sequences
        for sequence_idx, protein_sequence in enumerate(protein_sequences):
            fasta_string += f">protein|name={string.ascii_uppercase[sequence_idx]}\n{protein_sequence}\n"

        # Add the ligand
        fasta_string += f">ligand|name=Z\n{row['LIGAND_SMILES']}\n"

        output_path = os.path.join(args.output_folder, f"{row['PDB_CCD_ID']}.fasta")
        with open(output_path, "w") as f:
            f.write(fasta_string)


def generate_boltz_input(args: argparse.Namespace):
    """Generate Boltz input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    for _, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_dict = {"version": 1, "sequences": []}

        protein_sequences = row["PROTEIN_SEQUENCE"].split("|")
        # Skip if there are more than 20 chains
        if len(protein_sequences) > 20:
            print(f"Skipping {row['PDB_CCD_ID']} because it has more than 20 chains")
            continue
        
        # Add the protein sequences
        for sequence_idx, protein_sequence in enumerate(protein_sequences):
            task_dict["sequences"].append({
                "protein": {
                    "id": string.ascii_uppercase[sequence_idx],
                    "sequence": protein_sequence
                }
            })
        
        # Add the ligand
        task_dict["sequences"].append({
            "ligand": {
                "id": "Z",
                "smiles": row["LIGAND_SMILES"]
            }
        })

        # Save the task dict
        output_path = os.path.join(args.output_folder, f"{row['PDB_CCD_ID']}.yaml")
        with open(output_path, "w") as f:
            yaml.dump(task_dict, f)


def generate_rfaa_input(args: argparse.Namespace):
    """Generate RoseTTAFold-All-Atom input for a given docking data."""
    
    docking_data = pd.read_csv(args.input_file)
    for _, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_dict = {"defaults": ["base"], 
                     "job_name": row["PDB_CCD_ID"], 
                     "output_path": f"./predictions/{row['PDB_CCD_ID']}",
                     "protein_inputs": {},
                     "sm_inputs": {}}

        protein_sequences = row["PROTEIN_SEQUENCE"].split("|")
        # Skip if there are more than 20 chains
        if len(protein_sequences) > 20:
            print(f"Skipping {row['PDB_CCD_ID']} because it has more than 20 chains")
            continue

        # Add the protein sequences
        for sequence_idx, protein_sequence in enumerate(protein_sequences):
            chain_id = string.ascii_uppercase[sequence_idx]
            fasta_string = f">{row['PDB_CCD_ID']}|{chain_id}\n{protein_sequence}\n"
            fasta_path = os.path.join(args.output_folder, f"{row['PDB_CCD_ID']}_{chain_id}.fasta")
            with open(fasta_path, "w") as f:
                f.write(fasta_string)
            task_dict["protein_inputs"][chain_id] = {"fasta_file": os.path.abspath(fasta_path)}

        # Add the ligand
        task_dict["sm_inputs"]["Z"] = {"input": os.path.abspath(row["LIGAND_SDF_PATH"]), "input_type": "sdf"}

        # Save the task dict
        output_path = os.path.join(args.output_folder, f"{row['PDB_CCD_ID']}.yaml")
        with open(output_path, "w") as f:
            yaml.dump(task_dict, f)
    

def generate_dynamicbind_input(args: argparse.Namespace):
    """Generate DynaimcBind input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    for _, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        ligand_path = os.path.join(args.output_folder, f"{row['PDB_CCD_ID']}.csv")
        protein_path = os.path.join(args.output_folder, f"{row['PDB_CCD_ID']}.pdb")
        with open(ligand_path, "w") as f:
            f.write("ligand\n")
            f.write(f"{row['LIGAND_SMILES']}\n")
        shutil.copy(row['PROTEIN_PDB_PATH'], protein_path)


def generate_tankbind_input(args: argparse.Namespace):
    """Generate TankBind input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    input_data = docking_data[["PDB_CCD_ID", "PROTEIN_PDB_PATH", "LIGAND_SDF_PATH"]]
    input_data.to_csv(f"{args.output_folder}/data.csv", index=False)


def generate_diffdock_input(args: argparse.Namespace):
    """Generate DiffD input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    input_data = docking_data[["PDB_CCD_ID", "PROTEIN_PDB_PATH", "PROTEIN_SEQUENCE", "LIGAND_SMILES"]].copy()
    input_data.columns = ["complex_name", "protein_path", "protein_sequence", "ligand_description"]
    input_data["protein_path"] = input_data["protein_path"].apply(os.path.abspath)
    input_data.to_csv(f"{args.output_folder}/data.csv", index=False)


def generate_fabind_input(args: argparse.Namespace):
    """Generate FABind input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)    
    input_data = docking_data[["LIGAND_SMILES", "PDB_CCD_ID"]]
    input_data.to_csv(f"{args.output_folder}/ligand.csv", index=False)
    protein_dir = os.path.join(args.output_folder, "protein")
    os.makedirs(protein_dir, exist_ok=True)
    for _, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        protein_path = os.path.join(protein_dir, f"{row['PDB_CCD_ID']}.pdb")
        shutil.copy(row['PROTEIN_PDB_PATH'], protein_path)


def main(args: argparse.Namespace):
    os.makedirs(args.output_folder, exist_ok=True)
    print(f"Saving {args.model_type} input to {args.output_folder}")

    if args.model_type == "alphafold3":
        generate_alphafold3_input(args)
    elif args.model_type == "chai":
        generate_chai_input(args)
    elif args.model_type == "boltz":
        generate_boltz_input(args)
    elif args.model_type == "rfaa":
        generate_rfaa_input(args)
    elif args.model_type == "dynamicbind":
        generate_dynamicbind_input(args)
    elif args.model_type == "tankbind":
        generate_tankbind_input(args)
    elif args.model_type in ["diffdock", "diffdock_l"]:
        generate_diffdock_input(args)
    elif args.model_type == "fabind":
        generate_fabind_input(args)
    else:
        raise ValueError(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the alphafold3 input folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()

    main(args)
