# --------------------------------------------------------------------
# Following code adapted from (https://github.com/bytedance/Protenix)
# --------------------------------------------------------------------

import functools
import logging
from collections import defaultdict
from typing import Any, Optional, Union

import biotite
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import numpy as np
from biotite.structure import AtomArray

logger = logging.getLogger(__name__)

COMPONENTS_FILE = None


@functools.lru_cache
def biotite_load_ccd_cif() -> pdbx.CIFFile:
    """biotite load CCD components file

    Returns:
        pdbx.CIFFile: ccd components file
    """
    return pdbx.CIFFile.read(COMPONENTS_FILE)


def _map_central_to_leaving_groups(component) -> Optional[dict[str, list[list[str]]]]:
    """map each central atom (bonded atom) index to leaving atom groups in component (atom_array).

    Returns:
        dict[str, list[list[str]]]: central atom name to leaving atom groups (atom names).
    """
    comp = component.copy()
    # Eg: ions
    if comp.bonds is None:
        return {}
    central_to_leaving_groups = defaultdict(list)
    for c_idx in np.flatnonzero(~comp.leaving_atom_flag):
        bonds, _ = comp.bonds.get_bonds(c_idx)
        for l_idx in bonds:
            if comp.leaving_atom_flag[l_idx]:
                comp.bonds.remove_bond(c_idx, l_idx)
                group_idx = struc.find_connected(comp.bonds, l_idx)
                if not np.all(comp.leaving_atom_flag[group_idx]):
                    return None
                central_to_leaving_groups[comp.atom_name[c_idx]].append(
                    comp.atom_name[group_idx].tolist()
                )
    return central_to_leaving_groups


@functools.lru_cache
def get_component_atom_array(
    ccd_code: str, keep_leaving_atoms: bool = False, keep_hydrogens=False
) -> AtomArray:
    """get component atom array

    Args:
        ccd_code (str): ccd code
        keep_leaving_atoms (bool, optional): keep leaving atoms. Defaults to False.
        keep_hydrogens (bool, optional): keep hydrogens. Defaults to False.

    Returns:
        AtomArray: Biotite AtomArray of CCD component
            with additional attribute: leaving_atom_flag (bool)
    """
    ccd_cif = biotite_load_ccd_cif()
    if ccd_code not in ccd_cif:
        logger.warning(f"Warning: get_component_atom_array() can not parse {ccd_code}")
        return None
    try:
        comp = pdbx.get_component(ccd_cif, data_block=ccd_code, use_ideal_coord=True)
    except biotite.InvalidFileError as e:
        # Eg: UNL without atom.
        logger.warning(
            f"Warning: get_component_atom_array() can not parse {ccd_code} for {e}"
        )
        return None
    atom_category = ccd_cif[ccd_code]["chem_comp_atom"]
    leaving_atom_flag = atom_category["pdbx_leaving_atom_flag"].as_array()
    comp.set_annotation("leaving_atom_flag", leaving_atom_flag == "Y")

    for atom_id in ["alt_atom_id", "pdbx_component_atom_id"]:
        comp.set_annotation(atom_id, atom_category[atom_id].as_array())
    if not keep_leaving_atoms:
        comp = comp[~comp.leaving_atom_flag]
    if not keep_hydrogens:
        # EG: ND4
        comp = comp[~np.isin(comp.element, ["H", "D"])]

    # Map central atom index to leaving group (atom_indices) in component (atom_array).
    comp.central_to_leaving_groups = _map_central_to_leaving_groups(comp)
    if comp.central_to_leaving_groups is None:
        logger.warning(
            f"Warning: ccd {ccd_code} has leaving atom group bond to more than one central atom, central_to_leaving_groups is None."
        )
    return comp


@functools.lru_cache(maxsize=None)
def get_one_letter_code(ccd_code: str) -> Union[str, None]:
    """get one_letter_code from CCD components file.

    normal return is one letter: ALA --> A, DT --> T
    unknown protein: X
    unknown DNA or RNA: N
    other unknown: None
    some ccd_code will return more than one letter:
    eg: XXY --> THG

    Args:
        ccd_code (str): _description_

    Returns:
        str: one letter code
    """
    ccd_cif = biotite_load_ccd_cif()
    if ccd_code not in ccd_cif:
        return None
    one = ccd_cif[ccd_code]["chem_comp"]["one_letter_code"].as_item()
    if one == "?":
        return None
    else:
        return one


@functools.lru_cache(maxsize=None)
def get_mol_type(ccd_code: str) -> str:
    """get mol_type from CCD components file.

    based on _chem_comp.type
    http://mmcif.rcsb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.type.html

    not use _chem_comp.pdbx_type, because it is not consistent with _chem_comp.type
    e.g. ccd 000 --> _chem_comp.type="NON-POLYMER" _chem_comp.pdbx_type="ATOMP"
    https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v5_next.dic/Items/_struct_asym.pdbx_type.html

    Args:
        ccd_code (str): ccd code

    Returns:
        str: mol_type, one of {"protein", "rna", "dna", "ligand"}
    """
    ccd_cif = biotite_load_ccd_cif()
    if ccd_code not in ccd_cif:
        return "ligand"

    link_type = ccd_cif[ccd_code]["chem_comp"]["type"].as_item().upper()

    if "PEPTIDE" in link_type and link_type != "PEPTIDE-LIKE":
        return "protein"
    if "DNA" in link_type:
        return "dna"
    if "RNA" in link_type:
        return "rna"
    return "ligand"


def get_all_ccd_code() -> list:
    """get all ccd code from components file"""
    ccd_cif = biotite_load_ccd_cif()
    return list(ccd_cif.keys())

# Modified from biotite to use consistent ccd components file
def _connect_inter_residue(
    atoms: AtomArray, residue_starts: np.ndarray
) -> struc.BondList:
    """
    Create a :class:`BondList` containing the bonds between adjacent
    amino acid or nucleotide residues.

    Parameters
    ----------
    atoms : AtomArray or AtomArrayStack
        The structure to create the :class:`BondList` for.
    residue_starts : ndarray, dtype=int
        Return value of
        ``get_residue_starts(atoms, add_exclusive_stop=True)``.

    Returns
    -------
    BondList
        A bond list containing all inter residue bonds.
    """

    bonds = []

    atom_names = atoms.atom_name
    res_names = atoms.res_name
    res_ids = atoms.res_id
    chain_ids = atoms.chain_id

    # Iterate over all starts excluding:
    #   - the last residue and
    #   - exclusive end index of 'atoms'
    for i in range(len(residue_starts) - 2):
        curr_start_i = residue_starts[i]
        next_start_i = residue_starts[i + 1]
        after_next_start_i = residue_starts[i + 2]

        # Check if the current and next residue is in the same chain
        if chain_ids[next_start_i] != chain_ids[curr_start_i]:
            continue
        # Check if the current and next residue
        # have consecutive residue IDs
        # (Same residue ID is also possible if insertion code is used)
        if res_ids[next_start_i] - res_ids[curr_start_i] > 1:
            continue

        # Get link type for this residue from RCSB components.cif
        curr_link = get_mol_type(res_names[curr_start_i])
        next_link = get_mol_type(res_names[next_start_i])

        if curr_link == "protein" and next_link in "protein":
            curr_connect_atom_name = "C"
            next_connect_atom_name = "N"
        elif curr_link in ["dna", "rna"] and next_link in ["dna", "rna"]:
            curr_connect_atom_name = "O3'"
            next_connect_atom_name = "P"
        else:
            # Create no bond if the connection types of consecutive
            # residues are not compatible
            continue

        # Index in atom array for atom name in current residue
        # Addition of 'curr_start_i' is necessary, as only a slice of
        # 'atom_names' is taken, beginning at 'curr_start_i'
        curr_connect_indices = np.where(
            atom_names[curr_start_i:next_start_i] == curr_connect_atom_name
        )[0]
        curr_connect_indices += curr_start_i

        # Index in atom array for atom name in next residue
        next_connect_indices = np.where(
            atom_names[next_start_i:after_next_start_i] == next_connect_atom_name
        )[0]
        next_connect_indices += next_start_i

        if len(curr_connect_indices) == 0 or len(next_connect_indices) == 0:
            # The connector atoms are not found in the adjacent residues
            # -> skip this bond
            continue

        bonds.append(
            (curr_connect_indices[0], next_connect_indices[0], struc.BondType.SINGLE)
        )

    return struc.BondList(atoms.array_length(), np.array(bonds, dtype=np.uint32))


def add_inter_residue_bonds(
    atom_array: AtomArray,
    exclude_struct_conn_pairs: bool = False,
    remove_far_inter_chain_pairs: bool = False,
) -> AtomArray:
    """
    add polymer bonds (C-N or O3'-P) between adjacent residues based on auth_seq_id.

    exclude_struct_conn_pairs: if True, do not add bond between adjacent residues already has non-standard polymer bonds
                  on atom C or N or O3' or P.

    remove_far_inter_chain_pairs: if True, remove inter chain (based on label_asym_id) bonds that are far away from each other.

    returns:
        AtomArray: Biotite AtomArray merged inter residue bonds into atom_array.bonds
    """
    res_starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
    inter_bonds = _connect_inter_residue(atom_array, res_starts)

    if atom_array.bonds is None:
        atom_array.bonds = inter_bonds
        return atom_array

    select_mask = np.ones(len(inter_bonds._bonds), dtype=bool)
    if exclude_struct_conn_pairs:
        for b_idx, (atom_i, atom_j, b_type) in enumerate(inter_bonds._bonds):
            atom_k = atom_i if atom_array.atom_name[atom_i] in ("N", "O3'") else atom_j
            bonds, types = atom_array.bonds.get_bonds(atom_k)
            if len(bonds) == 0:
                continue
            for b in bonds:
                if (
                    # adjacent residues
                    abs((res_starts <= b).sum() - (res_starts <= atom_k).sum()) == 1
                    and atom_array.chain_id[b] == atom_array.chain_id[atom_k]
                    and atom_array.atom_name[b] not in ("C", "P")
                ):
                    select_mask[b_idx] = False
                    break

    if remove_far_inter_chain_pairs:
        if not hasattr(atom_array, "label_asym_id"):
            logging.warning(
                "label_asym_id not found, far inter chain bonds will not be removed"
            )
        for b_idx, (atom_i, atom_j, b_type) in enumerate(inter_bonds._bonds):
            if atom_array.label_asym_id[atom_i] != atom_array.label_asym_id[atom_j]:
                coord_i = atom_array.coord[atom_i]
                coord_j = atom_array.coord[atom_j]
                if np.linalg.norm(coord_i - coord_j) > 2.5:
                    select_mask[b_idx] = False

    # filter out removed_inter_bonds from atom_array.bonds
    remove_bonds = inter_bonds._bonds[~select_mask]
    remove_mask = np.isin(atom_array.bonds._bonds[:, 0], remove_bonds[:, 0]) & np.isin(
        atom_array.bonds._bonds[:, 1], remove_bonds[:, 1]
    )
    atom_array.bonds._bonds = atom_array.bonds._bonds[~remove_mask]

    # merged normal inter_bonds into atom_array.bonds
    inter_bonds._bonds = inter_bonds._bonds[select_mask]
    atom_array.bonds = atom_array.bonds.merge(inter_bonds)
    return atom_array


def res_names_to_sequence(res_names: list[str]) -> str:
    """convert res_names to sequences {chain_id: canonical_sequence} based on CCD

    Return
        str: canonical_sequence
    """
    seq = ""
    for res_name in res_names:
        one = get_one_letter_code(res_name)
        one = "X" if one is None else one
        one = "X" if len(one) > 1 else one
        seq += one
    return seq
