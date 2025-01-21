# --------------------------------------------------------------------
# Following code adapted from (https://github.com/bytedance/Protenix)
# --------------------------------------------------------------------

import os
import string
import copy
import gzip
import logging
import argparse
import functools
import numpy as np
import pandas as pd
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from posex import ccd
from pathlib import Path
from itertools import product
from collections import Counter
from typing import Optional, Union
from biotite.structure.io.pdbx import convert as pdbx_convert
from biotite.structure import AtomArray, get_chain_starts, get_residue_starts


CRYSTALLIZATION_AIDS = (
    "SO4",
    "GOL",
    "EDO",
    "PO4",
    "ACT",
    "PEG",
    "DMS",
    "TRS",
    "PGE",
    "PG4",
    "FMT",
    "EPE",
    "MPD",
    "MES",
    "CD",
    "IOD",
)


# Standard residues (AlphaFold3 SI Talbe 13)
PRO_STD_RESIDUES = {
    "ALA": 0,
    "ARG": 1,
    "ASN": 2,
    "ASP": 3,
    "CYS": 4,
    "GLN": 5,
    "GLU": 6,
    "GLY": 7,
    "HIS": 8,
    "ILE": 9,
    "LEU": 10,
    "LYS": 11,
    "MET": 12,
    "PHE": 13,
    "PRO": 14,
    "SER": 15,
    "THR": 16,
    "TRP": 17,
    "TYR": 18,
    "VAL": 19,
    "UNK": 20,
}

RNA_STD_RESIDUES = {
    "A": 21,
    "G": 22,
    "C": 23,
    "U": 24,
    "N": 29,
}

DNA_STD_RESIDUES = {
    "DA": 25,
    "DG": 26,
    "DC": 27,
    "DT": 28,
    "DN": 30,
}

STD_RESIDUES = PRO_STD_RESIDUES | RNA_STD_RESIDUES | DNA_STD_RESIDUES


class MMCIFParser:
    def __init__(self, mmcif_file: Union[str, Path]) -> None:
        self.cif = self._parse(mmcif_file=mmcif_file)

    def _parse(self, mmcif_file: Union[str, Path]) -> pdbx.CIFFile:
        mmcif_file = Path(mmcif_file)
        if mmcif_file.suffix == ".gz":
            with gzip.open(mmcif_file, "rt") as f:
                cif_file = pdbx.CIFFile.read(f)
        else:
            with open(mmcif_file, "rt") as f:
                cif_file = pdbx.CIFFile.read(f)
        return cif_file

    def get_category_table(self, name: str) -> Union[pd.DataFrame, None]:
        if name not in self.cif.block:
            return None
        category = self.cif.block[name]
        category_dict = {k: column.as_array() for k, column in category.items()}
        return pd.DataFrame(category_dict, dtype=str)

    @staticmethod
    def mse_to_met(atom_array: AtomArray) -> AtomArray:
        """
        Ref: AlphaFold3 SI chapter 2.1
        MSE residues are converted to MET residues.

        Args:
            atom_array (AtomArray): Biotite AtomArray object.

        Returns:
            AtomArray: Biotite AtomArray object after converted MSE to MET.
        """
        mse = atom_array.res_name == "MSE"
        se = mse & (atom_array.atom_name == "SE")
        atom_array.atom_name[se] = "SD"
        atom_array.element[se] = "S"
        atom_array.res_name[mse] = "MET"
        atom_array.hetero[mse] = False
        return atom_array

    @functools.cached_property
    def methods(self) -> list[str]:
        """the methods to get the structure

        most of the time, methods only has one method, such as 'X-RAY DIFFRACTION',
        but about 233 entries have multi methods, such as ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'].

        Allowed Values:
        https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_exptl.method.html

        Returns:
            list[str]: such as ['X-RAY DIFFRACTION'], ['ELECTRON MICROSCOPY'], ['SOLUTION NMR', 'THEORETICAL MODEL'],
                ['X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'], ['ELECTRON MICROSCOPY', 'SOLUTION NMR'], etc.
        """
        if "exptl" not in self.cif.block:
            return []
        else:
            methods = self.cif.block["exptl"]["method"]
            return methods.as_array()

    def get_poly_res_names(
        self, atom_array: Optional[AtomArray] = None
    ) -> dict[str, list[str]]:
        """get 3-letter residue names by combining mmcif._entity_poly_seq and atom_array

        if ref_atom_array is None: keep first altloc residue of the same res_id based in mmcif._entity_poly_seq
        if ref_atom_array is provided: keep same residue of ref_atom_array.

        Returns
            dict[str, list[str]]: label_entity_id --> [res_ids, res_names]
        """
        entity_res_names = {}
        if atom_array is not None:
            # build entity_id -> res_id -> res_name for input atom array
            res_starts = struc.get_residue_starts(atom_array, add_exclusive_stop=False)
            for start in res_starts:
                entity_id = atom_array.label_entity_id[start]
                res_id = atom_array.res_id[start]
                res_name = atom_array.res_name[start]
                if entity_id in entity_res_names:
                    entity_res_names[entity_id][res_id] = res_name
                else:
                    entity_res_names[entity_id] = {res_id: res_name}

        # build reference entity atom array, including missing residues
        entity_poly_seq = self.get_category_table("entity_poly_seq")
        if entity_poly_seq is None:
            return {}

        poly_res_names = {}
        for entity_id, poly_type in self.entity_poly_type.items():
            chain_mask = entity_poly_seq.entity_id == entity_id
            seq_mon_ids = entity_poly_seq.mon_id[chain_mask].to_numpy(dtype=str)

            # replace all MSE to MET in _entity_poly_seq.mon_id
            seq_mon_ids[seq_mon_ids == "MSE"] = "MET"

            seq_nums = entity_poly_seq.num[chain_mask].to_numpy(dtype=int)

            if np.unique(seq_nums).size == seq_nums.size:
                # no altloc residues
                poly_res_names[entity_id] = seq_mon_ids
                continue

            # filter altloc residues, eg: 181 ALA (altloc A); 181 GLY (altloc B)
            select_mask = np.zeros(len(seq_nums), dtype=bool)
            matching_res_id = seq_nums[0]
            for i, res_id in enumerate(seq_nums):
                if res_id != matching_res_id:
                    continue

                res_name_in_atom_array = entity_res_names.get(entity_id, {}).get(res_id)
                if res_name_in_atom_array is None:
                    # res_name is mssing in atom_array,
                    # keep first altloc residue of the same res_id
                    select_mask[i] = True
                else:
                    # keep match residue to atom_array
                    if res_name_in_atom_array == seq_mon_ids[i]:
                        select_mask[i] = True

                if select_mask[i]:
                    matching_res_id += 1

            seq_mon_ids = seq_mon_ids[select_mask]
            seq_nums = seq_nums[select_mask]
            assert len(seq_nums) == max(seq_nums)
            poly_res_names[entity_id] = seq_mon_ids
        return poly_res_names

    def get_sequences(self, atom_array=None) -> dict:
        """get sequence by combining mmcif._entity_poly_seq and atom_array

        if ref_atom_array is None: keep first altloc residue of the same res_id based in mmcif._entity_poly_seq
        if ref_atom_array is provided: keep same residue of atom_array.

        Return
            Dict{str:str}: label_entity_id --> canonical_sequence
        """
        sequences = {}
        for entity_id, res_names in self.get_poly_res_names(atom_array).items():
            seq = ccd.res_names_to_sequence(res_names)
            sequences[entity_id] = seq
        return sequences

    @functools.cached_property
    def entity_poly_type(self) -> dict[str, str]:
        """
        Ref: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_entity_poly.type.html
        Map entity_id to entity_poly_type.

        Allowed Value:
        · cyclic-pseudo-peptide
        · other
        · peptide nucleic acid
        · polydeoxyribonucleotide
        · polydeoxyribonucleotide/polyribonucleotide hybrid
        · polypeptide(D)
        · polypeptide(L)
        · polyribonucleotide

        Returns:
            Dict: a dict of label_entity_id --> entity_poly_type.
        """
        entity_poly = self.get_category_table("entity_poly")
        if entity_poly is None:
            return {}

        return {i: t for i, t in zip(entity_poly.entity_id, entity_poly.type)}

    def filter_altloc(self, atom_array: AtomArray, altloc: str = "first"):
        """
        altloc: "first", "A", "B", "global_largest", etc

        Filter first alternate coformation (altloc) of a given AtomArray.
        - normally first altloc_id is 'A'
        - but in one case, first altloc_id is '1' in 6uwi.cif

        biotite v0.41 can not handle diff res_name at same res_id.
        For example, in 2pxs.cif, there are two res_name (XYG|DYG) at res_id 63,
        need to keep the first XYG.
        """
        if altloc == "all":
            return atom_array

        altloc_id = altloc
        if altloc == "first":
            letter_altloc_ids = np.unique(atom_array.label_alt_id)
            if len(letter_altloc_ids) == 1 and letter_altloc_ids[0] == ".":
                return atom_array
            letter_altloc_ids = letter_altloc_ids[letter_altloc_ids != "."]
            altloc_id = np.sort(letter_altloc_ids)[0]

        return atom_array[np.isin(atom_array.label_alt_id, [altloc_id, "."])]

    @staticmethod
    def replace_auth_with_label(atom_array: AtomArray) -> AtomArray:
        # fix issue https://github.com/biotite-dev/biotite/issues/553
        atom_array.chain_id = atom_array.label_asym_id

        # reset ligand res_id
        res_id = copy.deepcopy(atom_array.label_seq_id)
        chain_starts = get_chain_starts(atom_array, add_exclusive_stop=True)
        for chain_start, chain_stop in zip(chain_starts[:-1], chain_starts[1:]):
            if atom_array.label_seq_id[chain_start] != ".":
                continue
            else:
                res_starts = get_residue_starts(
                    atom_array[chain_start:chain_stop], add_exclusive_stop=True
                )
                num = 1
                for res_start, res_stop in zip(res_starts[:-1], res_starts[1:]):
                    res_id[chain_start:chain_stop][res_start:res_stop] = num
                    num += 1

        atom_array.res_id = res_id.astype(int)
        return atom_array

    def get_structure(
        self,
        altloc: str = "first",
        model: int = 1,
        bond_lenth_threshold: Union[float, None] = 2.4,
    ) -> AtomArray:
        """
        Get an AtomArray created by bioassembly of MMCIF.

        altloc: "first", "all", "A", "B", etc
        model: the model number of the structure.
        bond_lenth_threshold: the threshold of bond length. If None, no filter will be applied.
                              Default is 2.4 Angstroms.

        Returns:
            AtomArray: Biotite AtomArray object created by bioassembly of MMCIF.
        """
        use_author_fields = True
        extra_fields = ["label_asym_id", "label_entity_id", "auth_asym_id"]  # chain
        extra_fields += ["label_seq_id", "auth_seq_id"]  # residue
        atom_site_fields = {
            "occupancy": "occupancy",
            "pdbx_formal_charge": "charge",
            "B_iso_or_equiv": "b_factor",
            "label_alt_id": "label_alt_id",
        }  # atom
        for atom_site_name, alt_name in atom_site_fields.items():
            if atom_site_name in self.cif.block["atom_site"]:
                extra_fields.append(alt_name)

        block = self.cif.block

        extra_fields = set(extra_fields)

        atom_site = block.get("atom_site")

        models = atom_site["pdbx_PDB_model_num"].as_array(np.int32)
        model_starts = pdbx_convert._get_model_starts(models)
        model_count = len(model_starts)

        if model == 0:
            raise ValueError("The model index must not be 0")
        # Negative models mean model indexing starting from last model

        model = model_count + model + 1 if model < 0 else model
        if model > model_count:
            raise ValueError(
                f"The file has {model_count} models, "
                f"the given model {model} does not exist"
            )

        model_atom_site = pdbx_convert._filter_model(atom_site, model_starts, model)
        # Any field of the category would work here to get the length
        model_length = model_atom_site.row_count
        atoms = AtomArray(model_length)

        atoms.coord[:, 0] = model_atom_site["Cartn_x"].as_array(np.float32)
        atoms.coord[:, 1] = model_atom_site["Cartn_y"].as_array(np.float32)
        atoms.coord[:, 2] = model_atom_site["Cartn_z"].as_array(np.float32)

        atoms.box = pdbx_convert._get_box(block)

        # The below part is the same for both, AtomArray and AtomArrayStack
        pdbx_convert._fill_annotations(
            atoms, model_atom_site, extra_fields, use_author_fields
        )

        bonds = struc.connect_via_residue_names(atoms, inter_residue=False)
        if "struct_conn" in block:
            conn_bonds = pdbx_convert._parse_inter_residue_bonds(
                model_atom_site, block["struct_conn"]
            )
            coord1 = atoms.coord[conn_bonds._bonds[:, 0]]
            coord2 = atoms.coord[conn_bonds._bonds[:, 1]]
            dist = np.linalg.norm(coord1 - coord2, axis=1)
            if bond_lenth_threshold is not None:
                conn_bonds._bonds = conn_bonds._bonds[dist < bond_lenth_threshold]
            bonds = bonds.merge(conn_bonds)
        atoms.bonds = bonds

        atom_array = self.filter_altloc(atoms, altloc=altloc)

        # inference inter residue bonds based on res_id (auth_seq_id) and label_asym_id.
        atom_array = ccd.add_inter_residue_bonds(
            atom_array,
            exclude_struct_conn_pairs=True,
            remove_far_inter_chain_pairs=True,
        )

        # use label_seq_id to match seq and structure
        atom_array = self.replace_auth_with_label(atom_array)

        # inference inter residue bonds based on new res_id (label_seq_id).
        # the auth_seq_id is not reliable, some are discontinuous (8bvh), some with insertion codes (6ydy).
        atom_array = ccd.add_inter_residue_bonds(
            atom_array, exclude_struct_conn_pairs=True
        )
        return atom_array

    def expand_assembly(
        self, structure: AtomArray, assembly_id: str = "1"
    ) -> AtomArray:
        """
        Expand the given assembly to all chains
        copy from biotite.structure.io.pdbx.get_assembly

        Args:
            structure (AtomArray): The AtomArray of the structure to expand.
            assembly_id (str, optional): The assembly ID in mmCIF file. Defaults to "1".
                                         If assembly_id is "all", all assemblies will be returned.

        Returns:
            AtomArray: The assembly AtomArray.
        """
        block = self.cif.block

        try:
            assembly_gen_category = block["pdbx_struct_assembly_gen"]
        except KeyError:
            logging.info(
                "File has no 'pdbx_struct_assembly_gen' category, return original structure."
            )
            return structure

        try:
            struct_oper_category = block["pdbx_struct_oper_list"]
        except KeyError:
            logging.info(
                "File has no 'pdbx_struct_oper_list' category, return original structure."
            )
            return structure

        assembly_ids = assembly_gen_category["assembly_id"].as_array(str)

        if assembly_id != "all":
            if assembly_id is None:
                assembly_id = assembly_ids[0]
            elif assembly_id not in assembly_ids:
                raise KeyError(f"File has no Assembly ID '{assembly_id}'")

        ### Calculate all possible transformations
        transformations = pdbx_convert._get_transformations(struct_oper_category)

        ### Get transformations and apply them to the affected asym IDs
        assembly = None
        assembly_1_mask = []
        for id, op_expr, asym_id_expr in zip(
            assembly_gen_category["assembly_id"].as_array(str),
            assembly_gen_category["oper_expression"].as_array(str),
            assembly_gen_category["asym_id_list"].as_array(str),
        ):
            # Find the operation expressions for given assembly ID
            # We already asserted that the ID is actually present
            if assembly_id == "all" or id == assembly_id:
                operations = pdbx_convert._parse_operation_expression(op_expr)
                asym_ids = asym_id_expr.split(",")
                # Filter affected asym IDs
                sub_structure = copy.deepcopy(
                    structure[..., np.isin(structure.label_asym_id, asym_ids)]
                )
                sub_assembly = pdbx_convert._apply_transformations(
                    sub_structure, transformations, operations
                )
                # Merge the chains with asym IDs for this operation
                # with chains from other operations
                if assembly is None:
                    assembly = sub_assembly
                else:
                    assembly += sub_assembly

                if id == "1":
                    assembly_1_mask.extend([True] * len(sub_assembly))
                else:
                    assembly_1_mask.extend([False] * len(sub_assembly))

        if assembly_id == "1" or assembly_id == "all":
            assembly.set_annotation("assembly_1", np.array(assembly_1_mask))
        return assembly


class AddAtomArrayAnnot(object):
    """
    The methods in this class are all designed to add annotations to an AtomArray
    without altering the information in the original AtomArray.
    """

    @staticmethod
    def add_token_mol_type(
        atom_array: AtomArray, sequences: dict[str, str]
    ) -> AtomArray:
        """
        Add molecule types in atom_arry.mol_type based on ccd pdbx_type.

        Args:
            atom_array (AtomArray): Biotite AtomArray object.
            sequences (dict[str, str]): A dict of label_entity_id --> canonical_sequence

        Return
            AtomArray: add atom_arry.mol_type = "protein" | "rna" | "dna" | "ligand"
        """
        mol_types = np.zeros(len(atom_array), dtype="U7")
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            entity_id = atom_array.label_entity_id[start]
            if entity_id not in sequences:
                # non-poly is ligand
                mol_types[start:stop] = "ligand"
                continue
            res_name = atom_array.res_name[start]

            mol_types[start:stop] = ccd.get_mol_type(res_name)

        atom_array.set_annotation("mol_type", mol_types)
        return atom_array

    @staticmethod
    def unique_chain_and_add_ids(atom_array: AtomArray) -> AtomArray:
        """
        Unique chain ID and add asym_id, entity_id, sym_id.
        Adds a number to the chain ID to make chain IDs in the assembly unique.
        Example: [A, B, A, B, C] ==> [A0, B0, A1, B1, C0]

        Args:
            atom_array (AtomArray): Biotite AtomArray object.

        Returns:
            AtomArray: Biotite AtomArray object with new annotations:
                - asym_id_int: np.array(int)
                - entity_id_int: np.array(int)
                - sym_id_int: np.array(int)
        """
        entity_id_uniq = np.sort(np.unique(atom_array.label_entity_id))
        entity_id_dict = {e: i for i, e in enumerate(entity_id_uniq)}
        asym_ids = np.zeros(len(atom_array), dtype=int)
        entity_ids = np.zeros(len(atom_array), dtype=int)
        sym_ids = np.zeros(len(atom_array), dtype=int)
        chain_ids = np.zeros(len(atom_array), dtype="U4")
        counter = Counter()
        start_indices = struc.get_chain_starts(atom_array, add_exclusive_stop=True)
        alphabet_list = list(string.ascii_uppercase)
        combinations = list(product(alphabet_list, repeat=2))
        combinations = sorted([''.join(pair) for pair in combinations])
        alphabet_list += combinations
        number_to_alpha_map = {n: a for n, a in enumerate(alphabet_list)}
        for i in range(len(start_indices) - 1):
            start_i = start_indices[i]
            stop_i = start_indices[i + 1]
            asym_ids[start_i:stop_i] = i

            entity_id = atom_array.label_entity_id[start_i]
            entity_ids[start_i:stop_i] = entity_id_dict[entity_id]

            sym_ids[start_i:stop_i] = counter[entity_id]
            counter[entity_id] += 1
            new_chain_id = f"{atom_array.chain_id[start_i]}{number_to_alpha_map[sym_ids[start_i]]}"
            chain_ids[start_i:stop_i] = new_chain_id

        atom_array.set_annotation("asym_id_int", asym_ids)
        atom_array.set_annotation("entity_id_int", entity_ids)
        atom_array.set_annotation("sym_id_int", sym_ids)
        atom_array.chain_id = chain_ids
        return atom_array


class Filter(object):
    """
    Ref: AlphaFold3 SI Chapter 2.5.4
    """

    @staticmethod
    def remove_hydrogens(atom_array: AtomArray) -> AtomArray:
        """remove hydrogens and deuteriums"""
        return atom_array[~np.isin(atom_array.element, ["H", "D"])]

    @staticmethod
    def remove_water(atom_array: AtomArray) -> AtomArray:
        """remove water (HOH) and deuterated water (DOD)"""
        return atom_array[~np.isin(atom_array.res_name, ["HOH", "DOD"])]

    @staticmethod
    def remove_element_X(atom_array: AtomArray) -> AtomArray:
        """
        remove element X
        following residues have element X:
        - UNX: unknown one atom or ion
        - UNL: unknown ligand, some atoms are marked as X
        - ASX: ASP/ASN ambiguous, two ambiguous atoms are marked as X, 6 entries in the PDB
        - GLX: GLU/GLN ambiguous, two ambiguous atoms are marked as X, 5 entries in the PDB
        """
        X_mask = np.zeros(len(atom_array), dtype=bool)
        starts = struc.get_residue_starts(atom_array, add_exclusive_stop=True)
        for start, stop in zip(starts[:-1], starts[1:]):
            res_name = atom_array.res_name[start]
            if res_name in ["UNX", "UNL"]:
                X_mask[start:stop] = True
        atom_array = atom_array[~X_mask]

        # map ASX to ASP, as ASP is more symmetric than ASN
        mask = atom_array.res_name == "ASX"
        atom_array.res_name[mask] = "ASP"
        atom_array.atom_name[mask & (atom_array.atom_name == "XD1")] = "OD1"
        atom_array.atom_name[mask & (atom_array.atom_name == "XD2")] = "OD2"
        atom_array.element[mask & (atom_array.element == "X")] = "O"

        # map GLX to GLU, as GLU is more symmetric than GLN
        mask = atom_array.res_name == "GLX"
        atom_array.res_name[mask] = "GLU"
        atom_array.atom_name[mask & (atom_array.atom_name == "XE1")] = "OE1"
        atom_array.atom_name[mask & (atom_array.atom_name == "XE2")] = "OE2"
        atom_array.element[mask & (atom_array.element == "X")] = "O"
        return atom_array

    @staticmethod
    def remove_crystallization_aids(
        atom_array: AtomArray, entity_poly_type: dict
    ) -> AtomArray:
        """remove crystallization aids, eg: SO4, GOL, etc.

        Only remove crystallization aids if the chain is not polymer.

        Ref: AlphaFold3 SI Chapter 2.5.4
        """
        non_aids_mask = ~np.isin(atom_array.res_name, CRYSTALLIZATION_AIDS)
        poly_mask = np.isin(atom_array.label_entity_id, list(entity_poly_type.keys()))
        return atom_array[poly_mask | non_aids_mask]


def _atom_array_to_input_json(
    atom_array: AtomArray,
    parser: MMCIFParser
) -> dict:
    """
    Convert a Biotite AtomArray to a dict that can be used as input to the model.

    Args:
        atom_array (AtomArray): Biotite Atom array.
        parser (MMCIFParser): Instantiated Protenix MMCIFParer.
        assembly_id (str, optional): Assembly ID. Defaults to None.

    Returns:
        dict: Protenix input json dict.
    """
    # get sequences after modified AtomArray
    entity_seq = parser.get_sequences(atom_array)
    # add unique chain id
    atom_array = AddAtomArrayAnnot.unique_chain_and_add_ids(atom_array)
    # find polymer modifications
    entity_id_to_mod_list = {}
    for entity_id, res_names in parser.get_poly_res_names(atom_array).items():
        modifications_list = []
        for idx, res_name in enumerate(res_names):
            if res_name not in STD_RESIDUES:
                position = idx + 1
                modifications_list.append([res_name, position])
        if modifications_list:
            entity_id_to_mod_list[entity_id] = modifications_list

    chain_starts = get_chain_starts(atom_array, add_exclusive_stop=False)
    chain_starts_atom_array = atom_array[chain_starts]
    json_dict = {
        "sequences": [],
    }
    unique_label_entity_id = np.unique(atom_array.label_entity_id)
    for label_entity_id in unique_label_entity_id:
        entity_dict = {}
        asym_chains = chain_starts_atom_array[
            chain_starts_atom_array.label_entity_id == label_entity_id
        ]
        entity_type = parser.entity_poly_type.get(label_entity_id, "ligand")
        if entity_type != "ligand":
            if entity_type == "polypeptide(L)":
                entity_type = "protein"
            else:
                # DNA, RNA, DNA/RNA hybrid, polypeptide(D), etc.
                continue
            chain_ids = list(asym_chains.chain_id)
            chain_ids = chain_ids[0] if len(chain_ids) == 1 else chain_ids
            entity_dict["id"] = chain_ids
            sequence = entity_seq.get(label_entity_id)
            entity_dict["sequence"] = sequence
            # add PTM info
            if label_entity_id in entity_id_to_mod_list:
                modifications = entity_id_to_mod_list[label_entity_id]
                entity_dict["modifications"] = [
                    {"ptmType": mod_ccd_code, "ptmPosition": position}
                    for mod_ccd_code, position in modifications
                ]
            json_dict["sequences"].append({entity_type: entity_dict})
        else:
            pass
        
    return json_dict


def cif_to_seq(
    pdbid: str,
    cif_dir: str,
    assembly_id: str = "1",
    altloc: str = "first",
) -> dict:
    """
    Convert mmcif file to sequence dict.

    Args:
        pdbid (str): PDBID
        cif_dir (str): folder containing pdb entries (cif format)
        assembly_id (str, optional): Assembly ID. Defaults to "1
        altloc (str, optional): Altloc selection. Defaults to "first"

    Returns:
        dict: sequence dict.
    """
    mmcif_file = os.path.join(cif_dir, f"{pdbid}.cif")
    parser = MMCIFParser(mmcif_file)
    atom_array = parser.get_structure(altloc, model=1, bond_lenth_threshold=None)

    # remove HOH from entities
    atom_array = Filter.remove_water(atom_array)
    atom_array = Filter.remove_hydrogens(atom_array)
    atom_array = parser.mse_to_met(atom_array)
    atom_array = Filter.remove_element_X(atom_array)

    # remove crystallization_aids
    if any(["DIFFRACTION" in m for m in parser.methods]):
        atom_array = Filter.remove_crystallization_aids(
            atom_array, parser.entity_poly_type
        )

    if assembly_id is not None:
        # expand created AtomArray by expand bioassembly
        atom_array = parser.expand_assembly(atom_array, assembly_id)
    seq_dict = _atom_array_to_input_json(
        atom_array,
        parser
    )
    return pdbid, seq_dict


if __name__ == "__main__":
    ccd.COMPONENTS_FILE = "/data/dataset/posex/ccd/components.cif"
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif_file", type=str, required=True, help="The cif file to parse")
    args = parser.parse_args()
    print(cif_to_seq(args.cif_file))
