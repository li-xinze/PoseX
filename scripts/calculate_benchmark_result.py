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

    for pdb_ccd_id in docking_data["PDB_CCD_ID"]:
        if args.model_type == "chai" and pdb_ccd_id in ["8FAV_4Y5"]:
            print(f"Skipping {pdb_ccd_id} for Chai")
            continue
        
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

    bust_data = pd.DataFrame(bust_dict)
    print("Number of Benchmark Data: ", len(bust_data))
    save_folder = os.path.dirname(args.input_file)

    # Calculate posebusters result
    buster = PoseBusters(config="redock", top_n=None)
    bust_results = buster.bust_table(bust_data, full_report=True)
    bust_results["PDB_CCD_ID"] = bust_dict["PDB_CCD_ID"]
    bust_results.to_csv(os.path.join(save_folder, f"posebusters_benchmark_result_{args.model_type}.csv"), index=False)

    # Calculate accuracy    
    accuracy = len(bust_results[bust_results["rmsd_≤_2å"] == True]) / len(bust_results)
    print(f"RMSD ≤ 2 Å: {accuracy * 100:.2f}%")

    # Calculate posebusters test result
    test_data = bust_results[POSEBUSTER_TEST_COLUMNS].copy()
    test_data.loc[:, "pb_valid"] = test_data.iloc[:, 1:].all(axis=1)
    valid_data = test_data[test_data["rmsd_≤_2å"] & test_data["pb_valid"]]
    print(f"RMSD ≤ 2 Å and PB Valid: {len(valid_data) / len(bust_results) * 100:.2f}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--dataset_folder", type=str, required=True, help="Path to the dataset folder")
    parser.add_argument("--model_output_folder", type=str, required=True, help="Path to the model output folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()

    main(args)
