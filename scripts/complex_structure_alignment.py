import argparse
import logging
import os
from typing import Optional, Any, Tuple, Union, Dict

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB.Model import Model
from rdkit import Chem
from rdkit.Geometry import Point3D
from scipy import spatial
from scipy import spatial as spa
from scipy.spatial.transform import Rotation
from scipy.optimize import Bounds, minimize
from tqdm.auto import tqdm


POSEBUSTERS_INPUT_FILE = "data/microcyto/posebusters/posebusters_benchmark.csv"
POSEBUSTERS_INPUT_FOLDER = "data/posebusters/posebusters_benchmark_set"
ALPHAFOLD3_OUTPUT_FOLDER = "data/microcyto/posebusters/alphafold3/output"


logging.basicConfig(format="[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def load_receptor(input_path: str) -> Model:
    """Load a receptor from a given PDB filepath.

    Args:
        input_path (str): Path to the PDB file containing the receptor

    Returns:
        Model: BioPython model containing the receptor
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("id", input_path)
    receptor = structure[0]
    return receptor


def load_ligand(input_path: str, remove_hs: bool = False, sanitize: bool = False) -> Chem.Mol:
    """Load a ligand from a given SDF filepath.

    Args:
        input_path (str): Path to the SDF file containing the ligand

    Returns:
        Chem.Mol: RDKit molecule containing the ligand
    """
    supplier = Chem.SDMolSupplier(input_path, sanitize=False, removeHs=False)
    mol = supplier[0]

    if sanitize:
        Chem.SanitizeMol(mol)
    if remove_hs:
        mol = Chem.RemoveHs(mol, sanitize=sanitize)
    return mol


def extract_complex_structure(receptor: Model, ligand: Optional[Chem.Mol], filter_hetero_residues: bool = False) -> Dict[str, Any]:
    """Extract the structure coordinates of a protein-ligand complex.

    Args:
        receptor (Model): BioPython model containing the receptor
        ligand (Optional[Chem.Mol]): RDKit molecule containing the ligand
        filter_hetero_residues (bool, optional): Whether to filter out hetero-residues (default: `False`).

    Returns:
        Dict[str, Any]: Dictionary containing the receptor, atom coordinates, CA coordinates, N coordinates, C coordinates, and valid chain IDs.
    """

    # Extract ligand coordinates if provided
    ligand_coords = ligand.GetConformer().GetPositions() if ligand else None

    total_atom_coords, total_ca_coords, total_n_coords, total_c_coords = [], [], [], []
    total_chain_residue_counts, min_chain_ligand_distances, valid_chain_ids = [], [], []

    # Extract coordinates for each chain
    for chain in receptor:
        chain_atom_coords, chain_ca_coords, chain_n_coords, chain_c_coords = [], [], [], []
        invalid_residue_ids = []
        residue_count = 0

        for residue in chain:
            # Skip if residue is a water molecule or a hetero-residue
            if residue.get_resname() == "HOH":
                invalid_residue_ids.append(residue.get_id())
                continue
            if filter_hetero_residues and len(residue.get_id()[0]) > 1:
                invalid_residue_ids.append(residue.get_id())
                continue
            
            # Extract CA, N, and C atom coordinates
            ca, n, c = None, None, None
            residue_atom_coords = []
            for atom in residue:
                if atom.get_name() == "CA":
                    ca = list(atom.get_vector())
                elif atom.get_name() == "N":
                    n = list(atom.get_vector())
                elif atom.get_name() == "C":
                    c = list(atom.get_vector())
                residue_atom_coords.append(list(atom.get_vector()))
                
            if ca is None or n is None or c is None:
                invalid_residue_ids.append(residue.get_id())
                continue
            
            # Append valid residue coordinates
            chain_ca_coords.append(ca)
            chain_n_coords.append(n)
            chain_c_coords.append(c)
            chain_atom_coords.append(np.array(residue_atom_coords))
            residue_count += 1
        
        # Remove invalid residues
        for residue_id in invalid_residue_ids:
            chain.detach_child(residue_id)

        # Compute minimum distance between ligand and chain
        if len(chain_atom_coords) > 0 and ligand_coords is not None:
            min_chain_ligand_distance = np.min(spatial.distance.cdist(ligand_coords, np.concatenate(chain_atom_coords, axis=0)))
        else:
            min_chain_ligand_distance = np.inf
        min_chain_ligand_distances.append(min_chain_ligand_distance)

        # Append chain coordinates
        total_chain_residue_counts.append(residue_count)
        total_atom_coords.append(chain_atom_coords)
        total_ca_coords.append(np.array(chain_ca_coords))
        total_n_coords.append(np.array(chain_n_coords))
        total_c_coords.append(np.array(chain_c_coords))

        # Append valid chain IDs
        if residue_count > 0:
            valid_chain_ids.append(chain.get_id())

    # Handle case where no valid chains are found
    min_chain_ligand_distances = np.array(min_chain_ligand_distances)
    if len(valid_chain_ids) == 0:
        valid_chain_ids.append(np.argmin(min_chain_ligand_distances))

    valid_atom_coords, valid_ca_coords, valid_n_coords, valid_c_coords = [], [], [], []
    valid_chain_residue_counts, invalid_chain_ids = [], []

    # Extract coordinates for valid chains
    for i, chain in enumerate(receptor):
        # Skip invalid chains
        if chain.get_id() not in valid_chain_ids:
            invalid_chain_ids.append(chain.get_id())
            continue
        
        # Append valid chain coordinates
        valid_atom_coords.append(total_atom_coords[i])
        valid_ca_coords.append(total_ca_coords[i])
        valid_n_coords.append(total_n_coords[i])
        valid_c_coords.append(total_c_coords[i])
        valid_chain_residue_counts.append(total_chain_residue_counts[i])

    # Concatenate valid chain coordinates
    total_atom_coords = [item for sublist in valid_atom_coords for item in sublist]
    total_ca_coords = np.concatenate(valid_ca_coords, axis=0)
    total_n_coords = np.concatenate(valid_n_coords, axis=0)
    total_c_coords = np.concatenate(valid_c_coords, axis=0)
    
    # Remove invalid chains
    for invalid_chain_id in invalid_chain_ids:
        receptor.detach_child(invalid_chain_id)
    
    # Check that the number of CA, N, and C atoms match
    assert len(total_ca_coords) == len(total_n_coords) == len(total_c_coords), "Number of Ca atoms does not match N and C atoms."

    return {
        "receptor": receptor,
        "atom_coords": total_atom_coords,
        "ca_coords": total_ca_coords,
        "n_coords": total_n_coords,
        "c_coords": total_c_coords,
    }


def align_predicted_structure(
    smoothing_factor: float, 
    reference_ca_coords: np.ndarray, 
    predicted_ca_coords: np.ndarray, 
    reference_ligand_coords: Optional[np.ndarray], 
    return_rotation: bool = False
) -> Union[Tuple[Rotation, np.ndarray, np.ndarray], float]:
    """Perform an alignment of apo and holo protein structures and ligand coordinates using an optimized smoothing factor.

    Args:
        smoothing_factor (float): Smoothing factor controlling the alignment
        reference_ca_coords (np.ndarray): Array of CA atom coordinates for a reference protein structure
        predicted_ca_coords (np.ndarray): Array of CA atom coordinates for a predicted protein structure
        reference_ligand_coords (Optional[np.ndarray]): Array of ligand coordinates from a reference dataset
        return_rotation (bool, optional): Whether to return the rotation matrix and centroids (default: `False`)

    Returns:
        Union[Tuple[Rotation, np.ndarray, np.ndarray], float]: 
            If return_rotation is `True`, returns a tuple containing rotation matrix (`Rotation`), 
            centroid of CA atoms for a reference protein (`np.ndarray`), 
            and centroid of CA atoms for a prediction (`np.ndarray`). 
            If return_rotation is `False`, returns the inverse root mean square error of reciprocal distances (`float`).
    """
    # Compute centroids
    if reference_ligand_coords is not None:
        reference_dists = spa.distance.cdist(reference_ca_coords, reference_ligand_coords)
        weights = np.exp(-1 * smoothing_factor * np.amin(reference_dists, axis=1))
        reference_ca_centroid = np.sum(np.expand_dims(weights, axis=1) * reference_ca_coords, axis=0) / np.sum(weights)
        predicted_ca_centroid = np.sum(np.expand_dims(weights, axis=1) * predicted_ca_coords, axis=0) / np.sum(weights)
    else:
        weights = None
        reference_ca_centroid = np.mean(reference_ca_coords, axis=0)
        predicted_ca_centroid = np.mean(predicted_ca_coords, axis=0)

    # Center coordinates
    centered_reference_ca_coords = reference_ca_coords - reference_ca_centroid
    centered_predicted_ca_coords = predicted_ca_coords - predicted_ca_centroid

    # Compute rotation matrix
    rotation, _ = spa.transform.Rotation.align_vectors(centered_reference_ca_coords, centered_predicted_ca_coords, weights)
    if return_rotation:  # Return rotation matrix and centroids
        return rotation, reference_ca_centroid, predicted_ca_centroid
    
    # Compute inverse root mean square error of reciprocal distances
    if reference_ligand_coords is not None:
        centered_reference_ligand_coords = reference_ligand_coords - reference_ca_centroid
        aligned_predicted_ca_coords = rotation.apply(centered_predicted_ca_coords)
        aligned_predicted_reference_dists = spa.distance.cdist(aligned_predicted_ca_coords, centered_reference_ligand_coords)
        inv_r_rmse = np.sqrt(np.mean(((1 / reference_dists) - (1 / aligned_predicted_reference_dists)) ** 2))
    else:
        inv_r_rmse = np.nan
    return inv_r_rmse


def run_structure_alignment(predicted_protein_pdb: str, predicted_ligand_sdf: Optional[str], 
                            reference_protein_pdb: str, reference_ligand_sdf: str, model_type: str):
    """Run structure alignment for a predicted protein-ligand complex.

    Args:
        predicted_protein_pdb (str): Path to the predicted protein structure in PDB format
        predicted_ligand_sdf (Optional[str]): Path to the predicted ligand structure in SDF format
        reference_protein_pdb (str): Path to the reference protein structure in PDB format
        reference_ligand_sdf (str): Path to the reference ligand structure in SDF format
        model_type (str): Model type
    """
    # Load receptor and ligand
    predicted_receptor = load_receptor(predicted_protein_pdb)
    reference_receptor = load_receptor(reference_protein_pdb)
    predicted_ligand = load_ligand(predicted_ligand_sdf, remove_hs=True, sanitize=False)
    reference_ligand = load_ligand(reference_ligand_sdf, remove_hs=True, sanitize=False)

    assert predicted_ligand is not None, "Predicted ligand after loading is None"
    assert reference_ligand is not None, "Reference ligand after loading is None"

    # Extract CA coordinates
    if model_type != "alphafold3":
        predicted_ca_coords = extract_complex_structure(predicted_receptor, reference_ligand, filter_hetero_residues=False)["ca_coords"]
        reference_ca_coords = extract_complex_structure(reference_receptor, reference_ligand, filter_hetero_residues=False)["ca_coords"]
    else:
        predicted_ca_coords = extract_complex_structure(predicted_receptor, reference_ligand, filter_hetero_residues=True)["ca_coords"]
        reference_ca_coords = extract_complex_structure(reference_receptor, reference_ligand, filter_hetero_residues=True)["ca_coords"]

    # Extract ligand coordinates
    predicted_ligand_conformer = predicted_ligand.GetConformer()
    predicted_ligand_coords = predicted_ligand.GetConformer().GetPositions()
    reference_ligand_coords = reference_ligand.GetConformer().GetPositions()

    # Optimize smoothing factor
    try:    
        smoothing_factor = minimize(
            align_predicted_structure,
            [0.1],
            bounds=Bounds([0.0], [1.0]),
            args=(reference_ca_coords, predicted_ca_coords, reference_ligand_coords),
            tol=1e-8,
        ).x
    except Exception as e:
        print(f"Error optimizing smoothing factor for {reference_protein_pdb}: {e}")
        return

    # Align predicted structure
    rotation, reference_ca_centroid, predicted_ca_centroid = align_predicted_structure(
        smoothing_factor,
        reference_ca_coords,
        predicted_ca_coords,
        reference_ligand_coords,
        return_rotation=True,
    )

    # Transform and save predicted protein
    predicted_protein = PandasPdb().read_pdb(predicted_protein_pdb)
    predicted_protein_raw = (
        predicted_protein.df["ATOM"][["x_coord", "y_coord", "z_coord"]]
        .to_numpy()
        .squeeze()
        .astype(np.float32)
    )
    predicted_protein_align = rotation.apply(predicted_protein_raw - predicted_ca_centroid) + reference_ca_centroid
    predicted_protein.df["ATOM"][["x_coord", "y_coord", "z_coord"]] = predicted_protein_align
    predicted_protein.to_pdb(
        path=predicted_protein_pdb.replace(".pdb", "_aligned.pdb"),
        records=["ATOM"],
        gz=False,
    )

    # Transform and save predicted ligand
    predicted_ligand_aligned = rotation.apply(predicted_ligand_coords - predicted_ca_centroid) + reference_ca_centroid
    for i in range(predicted_ligand.GetNumAtoms()):
        x, y, z = predicted_ligand_aligned[i]
        predicted_ligand_conformer.SetAtomPosition(i, Point3D(x, y, z))
    with Chem.SDWriter(predicted_ligand_sdf.replace(".sdf", f"_aligned.sdf")) as f:
        f.write(predicted_ligand)


def main(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    for _, row in tqdm(docking_data.iterrows(), total=len(docking_data), desc="Processing Alignment"):
        pdb_ccd_id = row["PDB_CCD_ID"]
        if args.model_type == "alphafold3":
            pdb_ccd_id = pdb_ccd_id.lower()

        if not os.path.exists(os.path.join(args.model_output_folder, f"{pdb_ccd_id}")):
            print(f"Directory {pdb_ccd_id} does not exist")
            continue

        predicted_protein_pdb = os.path.join(args.model_output_folder, f"{pdb_ccd_id}", f"{pdb_ccd_id}_model_protein.pdb")
        predicted_ligand_sdf = os.path.join(args.model_output_folder, f"{pdb_ccd_id}", f"{pdb_ccd_id}_model_ligand.sdf")
        if not os.path.exists(predicted_protein_pdb) or not os.path.exists(predicted_ligand_sdf):
            print(f"Predicted protein or ligand does not exist for {pdb_ccd_id}")
            continue

        reference_protein_pdb = os.path.join(args.dataset_folder, f"{pdb_ccd_id.upper()}", f"{pdb_ccd_id.upper()}_protein.pdb")
        reference_ligand_sdf = os.path.join(args.dataset_folder, f"{pdb_ccd_id.upper()}", f"{pdb_ccd_id.upper()}_ligand.sdf")

        run_structure_alignment(predicted_protein_pdb=predicted_protein_pdb,
                                predicted_ligand_sdf=predicted_ligand_sdf,
                                reference_protein_pdb=reference_protein_pdb,
                                reference_ligand_sdf=reference_ligand_sdf,
                                model_type=args.model_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the posebusters input file")
    parser.add_argument("--dataset_folder", type=str, required=True, help="Path to the dataset folder")
    parser.add_argument("--model_output_folder", type=str, required=True, help="Path to the model output folder")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    args = parser.parse_args()
    
    main(args)
