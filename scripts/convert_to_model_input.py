import argparse
import json
import os
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


def main(args: argparse.Namespace):
    os.makedirs(args.output_folder, exist_ok=True)
    print(f"Saving {args.model_type} input to {args.output_folder}")

    if args.model_type == "alphafold3":
        generate_alphafold3_input(args)
    elif args.model_type == "chai":
        generate_chai_input(args)
    elif args.model_type == "boltz":
        generate_boltz_input(args)
    else:
        raise ValueError(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the alphafold3 input folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()

    main(args)
