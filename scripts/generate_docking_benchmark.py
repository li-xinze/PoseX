import argparse
import json
import os

import pandas as pd
from loguru import logger


def generate_posex_benchmark(args: argparse.Namespace):
    """Generate the docking benchmark for the PoseX dataset

    Args:
        args (argparse.Namespace): The input arguments
    """
    # Get the docking dataset input folder
    if args.dataset == "astex":
        input_folder = os.path.join(args.input_folder, "astex_diverse_set")
    elif args.dataset == "posex_self_dock":
        input_folder = os.path.join(args.input_folder, "posex_self_docking_set")
    else:
        input_folder = os.path.join(args.input_folder, "posex_cross_docking_set")

    # Get the PDB_CCD_IDs
    pdb_ccd_ids = os.listdir(input_folder)
    logger.info(f"Number of PoseX {args.dataset} Data: {len(pdb_ccd_ids)}")
    
    molecule_smiles_list, protein_sequence_list, pdb_path_list, sdf_path_list = [], [], [], []
    for pdb_ccd_id in pdb_ccd_ids:
        with open(os.path.join(input_folder, pdb_ccd_id, f"{pdb_ccd_id}.json"), "r") as f:
            input_dict = json.load(f)
        
        # Get the protein sequences and the SMILES of the ligand
        chain_list = []
        for entity in input_dict["sequences"]:
            if "protein" in entity:
                chain_list.append(entity["protein"]["sequence"])
            elif "ligand" in entity:
                molecule_smiles_list.append(entity["ligand"]["smiles"])
        if len(chain_list) > 20:
            logger.warning(f"Warning: {pdb_ccd_id} has more than 20 chains.")
        if sum(len(chain) for chain in chain_list) > 2000:
            logger.warning(f"Warning: {pdb_ccd_id} has more than 2000 amino acids.")
        protein_sequence_list.append("|".join(chain_list))

        # Get the PDB file path and the SDF file path
        pdb_path_list.append(os.path.join(input_folder, pdb_ccd_id, f"{pdb_ccd_id}_protein.pdb"))
        sdf_path_list.append(os.path.join(input_folder, pdb_ccd_id, f"{pdb_ccd_id}_ligand_start_conf.sdf"))

    # Create the output data
    output_data = pd.DataFrame({
        "PDB_CCD_ID": pdb_ccd_ids, 
        "LIGAND_SMILES": molecule_smiles_list, 
        "PROTEIN_SEQUENCE": protein_sequence_list, 
        "PROTEIN_PDB_PATH": pdb_path_list, 
        "LIGAND_SDF_PATH": sdf_path_list
    })

    # Save the output data to a CSV file
    output_path = os.path.join(args.output_folder, f"{args.dataset}_benchmark.csv")
    output_data.to_csv(output_path, index=False)
    logger.info(f"Saved PoseX {args.dataset} Benchmark to {output_path}")


def main(args: argparse.Namespace):
    # Check if the dataset is valid
    if args.dataset not in ["astex", "posex_self_dock", "posex_cross_dock"]:
        raise ValueError(f"Unknown dataset: {args.dataset}")
    
    # Check if the output folder exists, if not create it
    if not os.path.exists(args.output_folder):
        logger.info(f"Output folder {args.output_folder} does not exist, creating it.")
        os.makedirs(args.output_folder)

    # Generate the docking benchmark
    generate_posex_benchmark(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder containing the PoseX dataset")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the output folder")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name (astex or posex_self_dock or posex_cross_dock)")
    args = parser.parse_args()

    main(args)
