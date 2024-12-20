import argparse
import os
from collections import defaultdict

import pandas as pd
from posebusters import PoseBusters


POSEBUSTER_TEST_COLUMNS = [
    # accuracy #
    "rmsd_≤_2å",
    # chemical validity and consistency #
    "mol_pred_loaded",
    "mol_true_loaded",
    "mol_cond_loaded",
    "sanitization",
    "molecular_formula",
    "molecular_bonds",
    "tetrahedral_chirality",
    "double_bond_stereochemistry",
    # intramolecular validity #
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "internal_energy",
    # intermolecular validity #
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
]


def main(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    bust_dict = defaultdict(list)

    total_samples = 0
    for pdb_ccd_id in docking_data["PDB_CCD_ID"]:
        mol_true = os.path.join(args.dataset_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_ligands.sdf")
        mol_cond = os.path.join(args.dataset_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_protein.pdb")

        if args.model_type == "alphafold3":
            pdb_ccd_id = pdb_ccd_id.lower()
        mol_pred = os.path.join(args.model_output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand_aligned.sdf")
        if not os.path.exists(mol_pred):
            print(f"File {mol_pred} does not exist")
            continue

        bust_dict["PDB_CCD_ID"].append(pdb_ccd_id.upper())
        bust_dict["mol_pred"].append(mol_pred)
        bust_dict["mol_true"].append(mol_true)
        bust_dict["mol_cond"].append(mol_cond)
        total_samples += 1

    bust_data = pd.DataFrame(bust_dict)
    print("Number of Benchmark Data: ", total_samples)
    save_folder = os.path.dirname(args.input_file)

    # Calculate posebusters result
    buster = PoseBusters(config="redock", top_n=None)
    buster.config["loading"]["mol_true"]["load_all"] = True
    bust_results = buster.bust_table(bust_data, full_report=True)
    bust_results["PDB_CCD_ID"] = bust_dict["PDB_CCD_ID"]
    bust_results.to_csv(os.path.join(save_folder, f"{args.dataset}_benchmark_result_{args.model_type}.csv"), index=False)

    # Calculate accuracy    
    accuracy = len(bust_results[bust_results["rmsd_≤_2å"] == True]) / total_samples
    print(f"RMSD ≤ 2 Å: {accuracy * 100:.2f}%")

    # Calculate posebusters test result
    test_data = bust_results[POSEBUSTER_TEST_COLUMNS].copy()
    test_data.loc[:, "pb_valid"] = test_data.iloc[:, 1:].all(axis=1)
    valid_data = test_data[test_data["rmsd_≤_2å"] & test_data["pb_valid"]]
    print(f"RMSD ≤ 2 Å and PB Valid: {len(valid_data) / total_samples * 100:.2f}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--dataset_folder", type=str, required=True, help="Path to the dataset folder")
    parser.add_argument("--model_output_folder", type=str, required=True, help="Path to the model output folder")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()

    main(args)
