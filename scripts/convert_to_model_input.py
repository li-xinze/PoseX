import argparse
import csv
import json
import os
import shutil
import string

import pandas as pd
import yaml
from tqdm.auto import tqdm
from rdkit import Chem
import numpy as np


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

        os.makedirs(args.output_folder, exist_ok=True)
        with open(os.path.join(args.output_folder, f"{task_name}.json"), "w") as f:
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


def generate_deepdock_input(args: argparse.Namespace):
    """Generate DeepDock input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        task_dir = os.path.join(args.output_folder, task_name)
        crystal_ligand_path = os.path.join(os.path.dirname(ligand_path), f"{task_name}_ligand.sdf")

        os.makedirs(task_dir, exist_ok=True)
        shutil.copy(crystal_ligand_path, task_dir)
        shutil.copy(protein_path, task_dir)
        shutil.copy(ligand_path, task_dir)


def generate_neuralplexer_input(args: argparse.Namespace):
    """Generate NeuralPlexer input for a given docking data."""

    def remove_insertion_code(protein_path, task_dir):
        with open(protein_path) as f:
            lines = f.readlines()
        content = ""
        for line in lines:
            if not line.startswith("ATOM") or not line[26].isalpha():
                content += line
        output_path = os.path.join(task_dir, os.path.basename(protein_path))
        with open(output_path, "w") as fw:
            fw.write(content)

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        task_dir = os.path.join(args.output_folder, task_name)
        os.makedirs(task_dir, exist_ok=True)
        remove_insertion_code(protein_path, task_dir)
        shutil.copy(ligand_path, task_dir)


def generate_gnina_input(args: argparse.Namespace):
    """Generate Gnina input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        crystal_ligand_path = os.path.join(os.path.dirname(ligand_path), f"{task_name}_ligand.sdf")
        task_dir = os.path.join(args.output_folder, task_name)
        os.makedirs(task_dir, exist_ok=True)
        shutil.copy(protein_path, task_dir)
        shutil.copy(ligand_path, task_dir)
        shutil.copy(crystal_ligand_path, task_dir)


def generate_unimol_input(args: argparse.Namespace):
    """Generate Unimol_docking v2 input for a given docking data."""

    def calculated_docking_grid_sdf(ligand_path, json_path, pdb_ccd_id, add_size=10):
        os.makedirs(json_path, exist_ok=True)
        output_grid = os.path.join(json_path, pdb_ccd_id + '.json')
        add_size = add_size
        mol = Chem.MolFromMolFile(ligand_path, sanitize=False)
        coords = mol.GetConformer(0).GetPositions().astype(np.float32)
        min_xyz = [min(coord[i] for coord in coords) for i in range(3)]
        max_xyz = [max(coord[i] for coord in coords) for i in range(3)]
        center = np.mean(coords, axis=0)
        size = [abs(max_xyz[i] - min_xyz[i]) for i in range(3)]
        center_x, center_y, center_z = center
        size_x, size_y, size_z = size
        size_x = size_x + add_size
        size_y = size_y + add_size
        size_z = size_z + add_size
        grid_info = {
            "center_x": float(center_x),
            "center_y": float(center_y),
            "center_z": float(center_z),
            "size_x": float(size_x),
            "size_y": float(size_y),
            "size_z": float(size_z)
        }
        with open(output_grid, 'w') as f:
            json.dump(grid_info, f, indent=4)

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        task_dir = os.path.join(args.output_folder, task_name)
        crystal_ligand_path = os.path.join(os.path.dirname(ligand_path), f"{task_name}_ligand.sdf")
        calculated_docking_grid_sdf(crystal_ligand_path, task_dir, task_name)
        shutil.copy(protein_path, task_dir)
        shutil.copy(ligand_path, task_dir)


def generate_interformer_input(args: argparse.Namespace):
    """Generate Interformer input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        crystal_ligand_path = os.path.join(os.path.dirname(ligand_path), f"{task_name}_ligand.sdf")
        task_dir = os.path.join(args.output_folder, task_name, "raw")
        os.makedirs(task_dir, exist_ok=True)
        shutil.copy(protein_path, task_dir)
        shutil.copy(ligand_path, task_dir)
        shutil.copy(crystal_ligand_path, task_dir)


def generate_equibind_input(args: argparse.Namespace):
    """Generate EquiBind input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        task_dir = os.path.join(args.output_folder, task_name)
        os.makedirs(task_dir, exist_ok=True)
        shutil.copy(protein_path, task_dir)
        shutil.copy(ligand_path, task_dir)


def generate_protenix_input(args: argparse.Namespace):
    """Generate Protenix input for a given docking data."""

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
                "proteinChain": {
                    "sequence": protein_sequence,
                    "count": 1,
                },
            }
            task_dict["sequences"].append(sequence_dict)

        # Add the ligand
        task_dict["sequences"].append({
            "ligand": {
                "ligand": row["LIGAND_SMILES"],
                "count": 1
            }
        })

        with open(os.path.join(args.output_folder, f"{task_name}.json"), "w") as f:
            json.dump([task_dict], f, indent=2)


def generate_surfdock_input(args: argparse.Namespace):
    """Generate SurfDock input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        crystal_ligand_path = os.path.join(os.path.dirname(ligand_path), f"{task_name}_ligand.sdf")
        task_dir = os.path.join(args.output_folder, task_name)
        os.makedirs(task_dir, exist_ok=True)
        protein_output_path = os.path.join(task_dir, f"{task_name}_protein_processed.pdb")
        shutil.copy(protein_path, protein_output_path)
        shutil.copy(crystal_ligand_path, task_dir)


def generate_diffdock_pocket_input(args: argparse.Namespace):
    """Generate DiffDock_Pocket input for a given docking data."""

    docking_data = pd.read_csv(args.input_file)
    head = ["complex_name", "experimental_protein", "ligand", "pocket_center_x", "pocket_center_y", "pocket_center_z",
            "flexible_sidechains"]
    data = [head]
    for data_idx, row in tqdm(docking_data.iterrows(), total=len(docking_data)):
        task_name = row["PDB_CCD_ID"]
        protein_path = row["PROTEIN_PDB_PATH"]
        ligand_path = row["LIGAND_SDF_PATH"]
        crystal_ligand_path = os.path.join(os.path.dirname(ligand_path), f"{task_name}_ligand.sdf")
        task_dir = os.path.join(args.output_folder, task_name)
        os.makedirs(task_dir, exist_ok=True)
        with open(protein_path) as f:
            lines = f.readlines()
        with open(os.path.join(task_dir, os.path.basename(protein_path)), "w") as fw:
            for line in lines:
                if not line.startswith("HETATM"):
                    fw.write(line)
        shutil.copy(crystal_ligand_path, task_dir)

        output_protein_path = os.path.join(task_dir, f"{task_name}_protein.pdb")
        output_ligand_path = os.path.join(task_dir, f"{task_name}_ligand.sdf")
        data.append(
            [task_name, os.path.abspath(output_protein_path), os.path.abspath(output_ligand_path), "", "", "", ""])
        with open(f"{args.output_folder}/example.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerows(data)


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
    elif args.model_type == "deepdock":
        generate_deepdock_input(args)
    elif args.model_type == "neuralplexer":
        generate_neuralplexer_input(args)
    elif args.model_type == "gnina":
        generate_gnina_input(args)
    elif args.model_type == "unimol":
        generate_unimol_input(args)
    elif args.model_type == "interformer":
        generate_interformer_input(args)
    elif args.model_type == "equibind":
        generate_equibind_input(args)
    elif args.model_type == "protenix":
        generate_protenix_input(args)
    elif args.model_type == "surfdock":
        generate_surfdock_input(args)
    elif args.model_type == "diffdock_pocket":
        generate_diffdock_pocket_input(args)
    else:
        raise ValueError(f"Unsupported model type: {args.model_type}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the alphafold3 input folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()

    main(args)
