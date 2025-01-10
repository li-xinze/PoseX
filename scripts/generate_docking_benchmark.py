import argparse
import os
from typing import List, Optional

import requests
import pandas as pd
from Bio.PDB import PDBParser
from rdkit import Chem
from tqdm.auto import tqdm

# RCSB PDB Entry URL
PDB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/"

# Normal amino acids residues
RAW_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# Specific amino acid residues or post-translational modifications
SPE_THREE_TO_ONE = {
    '2AS': 'D', '3AH': 'H', '5HP': 'E', 'ACL': 'R', 'AGM': 'R', 'AIB': 'A', 'ALM': 'A', 'ALO': 'T', 'ALY': 'K', 
    'ARM': 'R', 'ASA': 'D', 'ASB': 'D', 'ASK': 'D', 'ASL': 'D', 'ASQ': 'D', 'AYA': 'A', 'BCS': 'C', 'BHD': 'D', 
    'BMT': 'T', 'BNN': 'A', 'BUC': 'C', 'BUG': 'L', 'C5C': 'C', 'C6C': 'C', 'CAS': 'C', 'CCS': 'C', 'CEA': 'C', 
    'CGU': 'E', 'CHG': 'A', 'CLE': 'L', 'CME': 'C', 'CSD': 'A', 'CSO': 'C', 'CSP': 'C', 'CSS': 'C', 'CSW': 'C', 
    'CSX': 'C', 'CXM': 'M', 'CY1': 'C', 'CY3': 'C', 'CYG': 'C', 'CYM': 'C', 'CYQ': 'C', 'DAH': 'F', 'DAL': 'A', 
    'DAR': 'R', 'DAS': 'D', 'DCY': 'C', 'DGL': 'E', 'DGN': 'Q', 'DHA': 'A', 'DHI': 'H', 'DIL': 'I', 'DIV': 'V', 
    'DLE': 'L', 'DLY': 'K', 'DNP': 'A', 'DPN': 'F', 'DPR': 'P', 'DSN': 'S', 'DSP': 'D', 'DTH': 'T', 'DTR': 'W', 
    'DTY': 'Y', 'DVA': 'V', 'EFC': 'C', 'FLA': 'A', 'FME': 'M', 'GGL': 'E', 'GL3': 'G', 'GLZ': 'G', 'GMA': 'E', 
    'GSC': 'G', 'HAC': 'A', 'HAR': 'R', 'HIC': 'H', 'HIP': 'H', 'HMR': 'R', 'HPQ': 'F', 'HTR': 'W', 'HYP': 'P', 
    'IAS': 'D', 'IIL': 'I', 'IYR': 'Y', 'KCX': 'K', 'LLP': 'K', 'LLY': 'K', 'LTR': 'W', 'LYM': 'K', 'LYZ': 'K', 
    'MAA': 'A', 'MEN': 'N', 'MHS': 'H', 'MIS': 'S', 'MLE': 'L', 'MPQ': 'G', 'MSA': 'G', 'MSE': 'M', 'MVA': 'V', 
    'NEM': 'H', 'NEP': 'H', 'NLE': 'L', 'NLN': 'L', 'NLP': 'L', 'NMC': 'G', 'OAS': 'S', 'OCS': 'C', 'OMT': 'M', 
    'PAQ': 'Y', 'PCA': 'E', 'PEC': 'C', 'PHI': 'F', 'PHL': 'F', 'PR3': 'C', 'PRR': 'A', 'PTR': 'Y', 'PYX': 'C', 
    'SAC': 'S', 'SAR': 'G', 'SCH': 'C', 'SCS': 'C', 'SCY': 'C', 'SEL': 'S', 'SEP': 'S', 'SET': 'S', 'SHC': 'C', 
    'SHR': 'K', 'SMC': 'C', 'SOC': 'C', 'STY': 'Y', 'SVA': 'S', 'TIH': 'A', 'TPL': 'W', 'TPO': 'T', 'TPQ': 'A', 
    'TRG': 'K', 'TRO': 'W', 'TYB': 'Y', 'TYI': 'Y', 'TYQ': 'Y', 'TYS': 'Y', 'TYY': 'Y', 'I2M': 'I', 'MGN': 'Q',
    'MLY': 'K', '4OG': 'W', 'DYA': 'D'
}


def get_pdb_release_date(pdb_id: str) -> Optional[str]:
    """Get the release date of a PDB entry

    Args:
        pdb_id (str): PDB ID

    Returns:
        Optional[str]: Release date of the PDB entry
    """
    response = requests.get(PDB_ENTRY_URL + pdb_id)

    release_date = None
    if response.status_code == 200:
        data = response.json()
        release_date = data.get('rcsb_accession_info', {}).get('initial_release_date', None)
        if release_date:
            release_date = release_date.split("T")[0]

    return release_date


def get_molecule_smiles(input_file: str) -> str:
    """Read the smiles from a molecule file

    Args:
        input_file (str): Path to the molecule file

    Returns:
        str: SMILES of the molecule
    """
    mol = Chem.SDMolSupplier(input_file, sanitize=False, removeHs=False)[0]
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol, sanitize=True)
    return Chem.MolToSmiles(mol)


def get_protein_sequences(input_file: str) -> List[str]:
    """Get the protein sequences from a PDB file

    Args:
        input_file (str): Path to the PDB file

    Returns:
        List[str]: Protein sequences
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("id", input_file)[0]

    def check_residue(residue) -> bool:
        """Check if the residue is an amino acid or a post-translational modification"""
        if residue.get_resname() == "HOH":
            return False
        
        ca, n, c = None, None, None
        for atom in residue:
            if atom.get_name() == "CA":
                ca = atom.get_vector()
            elif atom.get_name() == "N":
                n = atom.get_vector()
            elif atom.get_name() == "C":
                c = atom.get_vector()
        if ca is None or n is None or c is None:
            return False
        
        return True

    protein_sequences = []
    for chain in structure:
        chain_sequence = ""
        
        for residue in chain:
            if not check_residue(residue):
                continue
            
            resname = residue.get_resname()
            # Normal amino acids
            if resname in RAW_THREE_TO_ONE:
                chain_sequence += RAW_THREE_TO_ONE[resname]
            # Specific amino acid residues or post-translational modifications
            elif resname in SPE_THREE_TO_ONE:
                chain_sequence += SPE_THREE_TO_ONE[resname]
                print(f"Warning: {input_file} {chain.id} {residue.id} {resname} is not a normal amino acid residue, replaced with {SPE_THREE_TO_ONE[resname]}.")
            # Unknown residues
            else:
                chain_sequence += "-"
                print(f"Warning: {input_file} {chain.id} {residue.id} {resname} cannot be recognized, replaced with a dash (i.e., `-`).")
        
        # Skip chains with less than 10 amino acids
        if len(chain_sequence) <= 10:
            print(f"Warning: {input_file} {chain.id} has less than 10 amino acids, skipping this chain.")
            continue

        protein_sequences.append(chain_sequence)
    
    return protein_sequences


def generate_posebusters_benchmark(args: argparse.Namespace):
    """Generate the Posebusters benchmark"""
    pdb_ccd_ids, pdb_ids, ccd_ids = [], [], []
    with open(os.path.join(args.input_folder, "posebusters_pdb_ccd_ids.txt"), "r") as f:
        for line in f:
            line = line.strip()
            pdb_ccd_ids.append(line)
            pdb_ids.append(line.split("_")[0])
            ccd_ids.append(line.split("_")[1])
    print("Number of Raw Posebusters Data:", len(pdb_ccd_ids))

    docking_data = pd.DataFrame({"PDB_CCD_ID": pdb_ccd_ids, "PDB_ID": pdb_ids, "CCD_ID": ccd_ids})

    pdb_ccd_ids, molecule_smiles_list, protein_sequence_list, pdb_path_list, sdf_path_list = [], [], [], [], []
    for pdb_ccd_id in docking_data["PDB_CCD_ID"]:
        data_folder = os.path.join(args.input_folder, f"posebusters_benchmark_set/{pdb_ccd_id}")

        # Get the protein sequences
        protein_file = os.path.join(data_folder, f"{pdb_ccd_id}_protein.pdb")
        protein_sequences = "|".join([seq for seq in get_protein_sequences(protein_file) if len(seq) > 0])
        if len(protein_sequences) > 1500:
            print(f"Warning: {pdb_ccd_id} has a protein sequence length of {len(protein_sequences)}, which is longer than 1500. Skipping this data.")
            continue
        if "-" in protein_sequences:
            print(f"Warning: {pdb_ccd_id} has a protein sequence containing a dash (i.e., `-`). Skipping this data.")
            continue

        # Get the SMILES of the ligand
        try:
            ligand_file = os.path.join(data_folder, f"{pdb_ccd_id}_ligand_start_conf.sdf")
            molecule_smiles = get_molecule_smiles(ligand_file)
        except Exception as e:
            print(f"Warning: {pdb_ccd_id} has an error when getting the SMILES of the ligand. Skipping this data.")
            continue
        
        # Append the data
        pdb_ccd_ids.append(pdb_ccd_id)
        molecule_smiles_list.append(molecule_smiles)
        protein_sequence_list.append(protein_sequences)
        pdb_path_list.append(protein_file)
        sdf_path_list.append(ligand_file)

    docking_data = docking_data[docking_data["PDB_CCD_ID"].isin(pdb_ccd_ids)].copy()
    docking_data["LIGAND_SMILES"] = molecule_smiles_list
    docking_data["PROTEIN_SEQUENCE"] = protein_sequence_list
    docking_data["PROTEIN_PDB_PATH"] = pdb_path_list
    docking_data["LIGAND_SDF_PATH"] = sdf_path_list
    print("Number of Filtered Posebusters Data:", len(docking_data))
 
    # Save the filtered data to a CSV file
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    docking_data.to_csv(os.path.join(args.output_folder, "posebusters_benchmark.csv"), index=False)
    print("Saved Posebusters Benchmark to ", os.path.join(args.output_folder, "posebusters_benchmark.csv"))


def generate_astex_benchmark(args: argparse.Namespace):
    """Generate the AsteX benchmark"""
    pdb_ccd_ids, pdb_ids, ccd_ids = [], [], []
    with open(os.path.join(args.input_folder, "astex_diverse_set_ids.txt"), "r") as f:
        for line in f:
            line = line.strip()
            pdb_ccd_ids.append(line)
            pdb_ids.append(line.split("_")[0])
            ccd_ids.append(line.split("_")[1])
    print("Number of Raw Astex Data:", len(pdb_ccd_ids))

    docking_data = pd.DataFrame({"PDB_CCD_ID": pdb_ccd_ids, "PDB_ID": pdb_ids, "CCD_ID": ccd_ids})

    pdb_ccd_ids, molecule_smiles_list, protein_sequence_list, pdb_path_list, sdf_path_list = [], [], [], [], []
    for pdb_ccd_id in docking_data["PDB_CCD_ID"]:
        data_folder = os.path.join(args.input_folder, f"astex_diverse_set/{pdb_ccd_id}")

        # Get the protein sequences
        protein_file = os.path.join(data_folder, f"{pdb_ccd_id}_protein.pdb")
        protein_sequences = "|".join([seq for seq in get_protein_sequences(protein_file) if len(seq) > 0])
        if len(protein_sequences) > 1500:
            print(f"Warning: {pdb_ccd_id} has a protein sequence length of {len(protein_sequences)}, which is longer than 2500. Skipping this data.")
            continue
        if "-" in protein_sequences:
            print(f"Warning: {pdb_ccd_id} has a protein sequence containing a dash (i.e., `-`). Skipping this data.")
            continue

        # Get the SMILES of the ligand
        try:
            ligand_file = os.path.join(data_folder, f"{pdb_ccd_id}_ligand_start_conf.sdf")
            molecule_smiles = get_molecule_smiles(ligand_file)
        except Exception as e:
            print(f"Warning: {pdb_ccd_id} has an error when getting the SMILES of the ligand. Skipping this data.")
            continue
        
        # Append the data
        pdb_ccd_ids.append(pdb_ccd_id)
        molecule_smiles_list.append(molecule_smiles)
        protein_sequence_list.append(protein_sequences)
        pdb_path_list.append(protein_file)
        sdf_path_list.append(ligand_file)

    docking_data = docking_data[docking_data["PDB_CCD_ID"].isin(pdb_ccd_ids)].copy()
    docking_data["LIGAND_SMILES"] = molecule_smiles_list
    docking_data["PROTEIN_SEQUENCE"] = protein_sequence_list
    docking_data["PROTEIN_PDB_PATH"] = pdb_path_list
    docking_data["LIGAND_SDF_PATH"] = sdf_path_list
    print("Number of Filtered Astex Data:", len(docking_data))

    # Save the filtered data to a CSV file
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    docking_data.to_csv(os.path.join(args.output_folder, "astex_benchmark.csv"), index=False)
    print("Saved Astex Benchmark to ", os.path.join(args.output_folder, "astex_benchmark.csv"))


def generate_posex_benchmark(args: argparse.Namespace, mode: str):
    """Generate the PoseX benchmark"""
    data_root = os.path.join(args.input_folder, f"posex_{mode}_set")
    pdb_ccd_ids = os.listdir(data_root)
    print(f"Number of PoseX {mode} Data:", len(pdb_ccd_ids))
    docking_data = pd.DataFrame({"PDB_CCD_ID": pdb_ccd_ids})
    pdb_ccd_ids, molecule_smiles_list, protein_sequence_list, pdb_path_list, sdf_path_list = [], [], [], [], []
    for pdb_ccd_id in docking_data["PDB_CCD_ID"]:
        data_folder = os.path.join(data_root, pdb_ccd_id)
        # Get the protein sequences
        protein_file = os.path.join(data_folder, f"{pdb_ccd_id}_protein.pdb")
        protein_sequences = "|".join([seq for seq in get_protein_sequences(protein_file) if len(seq) > 0])
        if len(protein_sequences) > 1500:
            print(f"Warning: {pdb_ccd_id} has a protein sequence length of {len(protein_sequences)}, which is longer than 2500. Skipping this data.")
            continue
        if "-" in protein_sequences:
            print(f"Warning: {pdb_ccd_id} has a protein sequence containing a dash (i.e., `-`). Skipping this data.")
            continue
        # Get the SMILES of the ligand
        try:
            ligand_file = os.path.join(data_folder, f"{pdb_ccd_id}_ligand_start_conf.sdf")
            molecule_smiles = get_molecule_smiles(ligand_file)
        except Exception as e:
            print(f"Warning: {pdb_ccd_id} has an error when getting the SMILES of the ligand. Skipping this data.")
            continue
        
        # Append the data
        pdb_ccd_ids.append(pdb_ccd_id)
        molecule_smiles_list.append(molecule_smiles)
        protein_sequence_list.append(protein_sequences)
        pdb_path_list.append(protein_file)
        sdf_path_list.append(ligand_file)

    docking_data = docking_data[docking_data["PDB_CCD_ID"].isin(pdb_ccd_ids)].copy()
    docking_data["LIGAND_SMILES"] = molecule_smiles_list
    docking_data["PROTEIN_SEQUENCE"] = protein_sequence_list
    docking_data["PROTEIN_PDB_PATH"] = pdb_path_list
    docking_data["LIGAND_SDF_PATH"] = sdf_path_list
    print("Number of Filtered Astex Data:", len(docking_data))

    # Save the filtered data to a CSV file
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    benchmark_file_name = f"posex_{mode}_benchmark.csv"
    docking_data.to_csv(os.path.join(args.output_folder, benchmark_file_name), index=False)
    print(f"Saved PoseX {mode} Benchmark to ", os.path.join(args.output_folder, benchmark_file_name))


def main(args: argparse.Namespace):
    if args.dataset == "posebusters":
        generate_posebusters_benchmark(args)
    elif args.dataset == "astex":
        generate_astex_benchmark(args)
    elif args.dataset == "posex_self_dock":
        generate_posex_benchmark(args, mode="self_dock")
    elif args.dataset == "posex_cross_dock":
        generate_posex_benchmark(args, mode="cross_dock")
    else:
        raise ValueError(f"Unknown dataset: {args.dataset}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder containing the Posebusters dataset")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the output folder")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset type")
    args = parser.parse_args()

    main(args)
