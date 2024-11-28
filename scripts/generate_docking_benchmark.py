import argparse
import os
import requests
from typing import List, Optional

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
    "0JO": "L", "4OG": "E", "AGM": "R", "CME": "C", "CSO": "C",
    "DGL": "E", "DYA": "Y", "FME": "M", "GL3": "Q", "KCX": "K",
    "LCS": "C", "LLP": "K", "MGN": "G", "MSE": "M", "MHS": "C",
    "OCS": "C", "PLS": "K", "PYL": "O", "SAM": "M", "SEP": "S",
    "SEC": "U", "SMC": "C", "TP7": "W", "TPO": "T", "ASX": "B",
    "GLX": "Z", "XAA": "X", "XLE": "J",
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
        
        protein_sequences.append(chain_sequence)
    
    return protein_sequences


def main(args: argparse.Namespace):
    pdb_ccd_ids, pdb_ids, ccd_ids = [], [], []
    with open(os.path.join(args.input_folder, "posebusters_pdb_ccd_ids.txt"), "r") as f:
        for line in f:
            line = line.strip()
            pdb_ccd_ids.append(line)
            pdb_ids.append(line.split("_")[0])
            ccd_ids.append(line.split("_")[1])
    print("Number of Raw Posebusters Data:", len(pdb_ccd_ids))

    # Get the release dates of the PDB entries
    release_dates = [get_pdb_release_date(pdb_id) for pdb_id in tqdm(pdb_ids, desc="Getting Release Dates")]
    docking_data = pd.DataFrame({"PDB_CCD_ID": pdb_ccd_ids, "PDB_ID": pdb_ids, "CCD_ID": ccd_ids, "RELEASE_DATE": release_dates})
    # Filter the data to only include PDB-CCD pairs with a release date after 2021-10-01
    docking_data = docking_data[docking_data["RELEASE_DATE"] >= "2021-10-01"]
    # Filter out sequences that are too long
    docking_data = docking_data[~docking_data["PDB_CCD_ID"].isin(["8F4J_PHO", "7SUC_COM"])]
    print("Number of Filtered Posebusters Data:", len(docking_data))

    molecule_smiles_list, protein_sequence_list = [], []
    for pdb_ccd_id in docking_data["PDB_CCD_ID"]:
        data_folder = os.path.join(args.input_folder, f"posebusters_benchmark_set/{pdb_ccd_id}")
        
        # Get the SMILES of the ligand
        molecule_smiles = get_molecule_smiles(os.path.join(data_folder, f"{pdb_ccd_id}_ligand.sdf"))
        molecule_smiles_list.append(molecule_smiles)

        # Get the protein sequences
        protein_file = os.path.join(data_folder, f"{pdb_ccd_id}_protein.pdb")
        protein_sequences = "|".join([seq for seq in get_protein_sequences(protein_file) if len(seq) > 0])
        protein_sequence_list.append(protein_sequences)

    docking_data["LIGAND_SMILES"] = molecule_smiles_list
    docking_data["PROTEIN_SEQUENCE"] = protein_sequence_list

    # Save the filtered data to a CSV file
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    docking_data.to_csv(os.path.join(args.output_folder, "posebusters_benchmark.csv"), index=False)
    print("Saved Posebusters Benchmark to ", os.path.join(args.output_folder, "posebusters_benchmark.csv"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder containing the Posebusters dataset")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the output folder")
    args = parser.parse_args()

    main(args)
