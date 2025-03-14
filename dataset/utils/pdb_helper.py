from pathlib import Path
from typing import *

import numpy as np
from lxml import etree
from pdbfixer.pdbfixer import PDBFixer
from rdkit import Chem

PROTEIN_RESIDUES = [
    "ALA",
    "ASN",
    "CYS",
    "GLU",
    "HIS",
    "LEU",
    "MET",
    "PRO",
    "THR",
    "TYR",
    "ARG",
    "ASP",
    "GLN",
    "GLY",
    "ILE",
    "LYS",
    "PHE",
    "SER",
    "TRP",
    "VAL",
]
METALS = ("NA", "K", "CA", "MG", "FE", "ZN", "CU", "MN", "CO", "NI")
WATERS = ("HOH", "DOD")


def format_4letter(atom_name: str):
    output = None
    if len(atom_name) == 4:
        output = atom_name
    elif len(atom_name) == 3:
        output = " " + atom_name
    elif len(atom_name) == 2:
        output = " " + atom_name + " "
    elif len(atom_name) == 1:
        output = " " + atom_name + "  "
    else:
        raise ValueError()

    return output


def load_amber_xml() -> Dict[str, int]:
    charge_table = {}
    hydrogen_table = {}
    connect_table = {}
    protein_xml_fn = Path(__file__).parent.parent / "amberlib/protein.ff14SB.xml"
    # 加载XML文件
    tree = etree.parse(protein_xml_fn)
    root = tree.getroot()
    # 遍历XML文件，获取每个元素的属性值
    # 遍历XML树
    for element in root:
        # print(element.tag, element.attrib)
        if element.tag == "Residues":
            for residue_element in element:
                res_name = residue_element.attrib["name"]
                charge_table[res_name] = {}
                connect_table[res_name] = {}
                hydrogen_table[res_name] = {}
                for element in residue_element:
                    if element.tag == "Atom":
                        atom_name = element.attrib["name"]
                        if atom_name in ["NH2", "NZ"]:
                            charge_table[res_name][atom_name] = 1
                        elif (
                            len(res_name) == 4
                            and res_name[0] == "N"
                            and atom_name == "N"
                        ):
                            # 不加冒的N端
                            charge_table[res_name][atom_name] = 1

                        elif res_name[-3:] == "HIP" and atom_name == "ND1":
                            charge_table[res_name][atom_name] = 1

                        elif atom_name in ["OXT", "OD2", "OE2"]:
                            charge_table[res_name][atom_name] = -1
                        else:
                            charge_table[res_name][atom_name] = 0

                    if element.tag == "Bond":
                        atom1 = element.attrib["atomName1"]
                        atom2 = element.attrib["atomName2"]
                        if set([atom1, atom2]) == set(["C", "O"]):
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2
                        elif res_name[-3:] == "ARG" and set([atom1, atom2]) == set(
                            ["CZ", "NH2"]
                        ):
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2

                        elif res_name[-3:] in ["HIS", "HIE", "HIP"] and set(
                            [atom1, atom2]
                        ) in [set(["CG", "CD2"]), set(["CE1", "ND1"])]:
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2
                        elif res_name[-3:] == "HID" and set([atom1, atom2]) in [
                            set(["CG", "CD2"]),
                            set(["CE1", "NE2"]),
                        ]:
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2

                        elif res_name[-3:] in ["ASP", "ASN", "ASH"] and set(
                            [atom1, atom2]
                        ) == set(["CG", "OD1"]):
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2

                        elif res_name[-3:] in ["GLU", "GLN", "GLH"] and set(
                            [atom1, atom2]
                        ) == set(["CD", "OE1"]):
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2

                        elif res_name[-3:] in ["PHE", "TYR"] and set(
                            [atom1, atom2]
                        ) in [
                            set(["CG", "CD1"]),
                            set(["CE1", "CZ"]),
                            set(["CE2", "CD2"]),
                        ]:
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2

                        elif res_name[-3:] == "TRP" and set([atom1, atom2]) in [
                            set(["CG", "CD1"]),
                            set(["CE2", "CD2"]),
                            set(["CH2", "CZ2"]),
                            set(["CE3", "CZ3"]),
                        ]:
                            connect_table[res_name][(atom1, atom2)] = 2
                            connect_table[res_name][(atom2, atom1)] = 2

                        else:
                            connect_table[res_name][(atom1, atom2)] = 1
                            connect_table[res_name][(atom2, atom1)] = 1

                        # if res_name == 'NME':
                        #     pass
                        if atom2[0] == "H" and atom1[0] != "H":
                            if atom1 not in hydrogen_table[res_name]:
                                hydrogen_table[res_name][atom1] = []
                            hydrogen_table[res_name][atom1].append(atom2)
                        elif atom1[0] == "H" and atom2[0] != "H":
                            if atom2 not in hydrogen_table[res_name]:
                                hydrogen_table[res_name][atom2] = []
                            hydrogen_table[res_name][atom2].append(atom1)

    hydrogen_table["HOH"] = {}
    hydrogen_table["HOH"]["O"] = ["H1", "H2"]

    charge_table["HOH"] = {}
    charge_table["HOH"]["O"] = 0
    charge_table["HOH"]["H1"] = 0
    charge_table["HOH"]["H2"] = 0
    for metal_ in METALS:
        charge_table[metal_] = {}

    charge_table["NA"]["NA"] = 1
    charge_table["K"]["K"] = 1
    charge_table["CA"]["CA"] = 2
    charge_table["MG"]["MG"] = 2
    charge_table["FE"]["FE"] = 2
    # charge_table["FE2"]["FE2"] = 2
    charge_table["ZN"]["ZN"] = 2
    charge_table["CU"]["CU"] = 2
    charge_table["MN"]["MN"] = 2
    charge_table["CO"]["CO"] = 3
    charge_table["NI"]["NI"] = 2
    # charge_table["CL"]["CL"] = -1

    connect_table["HOH"] = {}
    connect_table["HOH"][("O", "H1")] = 1
    connect_table["HOH"][("H1", "O")] = 1
    connect_table["HOH"][("O", "H2")] = 1
    connect_table["HOH"][("H2", "O")] = 1

    hydrogen_table["HIS"] = hydrogen_table["HID"]
    charge_table["HIS"] = charge_table["HID"]
    connect_table["HIS"] = connect_table["HID"]

    return connect_table, charge_table, hydrogen_table


def add_conformer_to_mol(new_mol: Chem.Mol, pos) -> Chem.Mol:
    conf = Chem.Conformer(new_mol.GetNumAtoms())
    for idx in range(len(pos)):
        conf.SetAtomPosition(idx, pos[idx])
    # conf.SetPositions(pos)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conf, assignId=True)
    return new_mol


def add_h_to_receptor_mol(rec_mol: Chem.Mol, hydrogen_table: Dict[str, List[str]]):
    new_mol = Chem.AddHs(rec_mol, addCoords=True)
    residue_mapper = {}
    for atom in new_mol.GetAtoms():
        atom: Chem.Atom
        if atom.GetAtomicNum() != 1:
            mi = atom.GetPDBResidueInfo()
            res_key = [mi.GetResidueName(), int(mi.GetResidueNumber()), mi.GetChainId()]
            atom_name = mi.GetName().strip()

            tmp_mi = Chem.AtomPDBResidueInfo()
            tmp_mi.SetResidueName(f"{res_key[0]:>3s}")
            tmp_mi.SetResidueNumber(res_key[1])
            tmp_mi.SetChainId(res_key[2])
            if mi.GetInsertionCode() == " ":
                tmp_mi.SetInsertionCode(" ")
            else:
                tmp_mi.SetInsertionCode(mi.GetInsertionCode())

            if tuple(res_key) not in residue_mapper:
                residue_mapper[tuple(res_key)] = []

            residue_mapper[tuple(res_key)].append(atom.GetIdx())

            h_atoms: List[Chem.Atom] = [
                nbr_atom
                for nbr_atom in atom.GetNeighbors()
                if nbr_atom.GetAtomicNum() == 1
            ]
            if len(h_atoms) > 0:
                std_h_names = hydrogen_table[res_key[0]][atom_name]
                assert len(h_atoms) == len(std_h_names)
                for h_atom, h_name in zip(h_atoms, list(sorted(std_h_names))):
                    tmp_mi.SetName(format_4letter(h_name))
                    h_atom.SetMonomerInfo(tmp_mi)
                    residue_mapper[tuple(res_key)].append(h_atom.GetIdx())

    resort_list = []
    [resort_list.extend(values) for _, values in residue_mapper.items()]
    assert len(resort_list) == new_mol.GetNumAtoms()
    new_mol_h = Chem.RenumberAtoms(new_mol, resort_list)

    return new_mol_h
