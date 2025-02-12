from rdkit import Chem
import openmm.app as app
from openmm import unit, Vec3
from collections import defaultdict

import numpy as np

CURRENTLY_ACCEPTABLE_PROTEIN_RESIDUES = {
    "ALA": "A",
    "ASN": "N",
    "CYS": "C",
    "GLU": "E",
    "HIS": "H",
    "LEU": "L",
    "MET": "M",
    "PRO": "P",
    "THR": "T",
    "TYR": "Y",
    "ARG": "R",
    "ASP": "D",
    "GLN": "Q",
    "GLY": "G",
    "ILE": "I",
    "LYS": "K",
    "PHE": "F",
    "SER": "S",
    "TRP": "W",
    "VAL": "V",
    "ACE": "B",
    "NME": "Z",
    "HOH": "O",
}

AA2SMILES = {
    "ACE": "CC(=O)",
    "NME": "NC",
    "ALA": "C[C@H](N)C=O",
    "CYS": "N[C@H](C=O)CS",
    "ASP": "N[C@H](C=O)CC(=O)[O-]",
    "GLU": "N[C@H](C=O)CCC(=O)[O-]",
    "PHE": "N[C@H](C=O)Cc1ccccc1",
    "GLY": "NCC=O",
    "HIS": [
        "N([H])C([H])(C=O)C([H])([H])C1=C([H])N([H])C([H])=N1",
        "N(C(C(C1N([H])C([H])=NC=1[H])([H])[H])(C=O)[H])[H]",
        "NC(CC1N=CNC=1)C=O",
    ],
    "ILE": "CC[C@H](C)[C@H](N)C=O",
    "LYS": "[NH3+]CCCC[C@H](N)C=O",
    "LEU": "CC(C)C[C@H](N)C=O",
    "MET": "CSCC[C@H](N)C=O",
    "ASN": "NC(=O)C[C@H](N)C=O",
    "PRO": "O=C[C@@H]1CCCN1",
    "GLN": "NC(=O)CC[C@H](N)C=O",
    "ARG": "NC(=[NH2+])NCCC[C@H](N)C=O",
    "SER": "N[C@H](C=O)CO",
    "THR": "C[C@@H](O)[C@H](N)C=O",
    "VAL": "CC(C)[C@H](N)C=O",
    "TRP": "N[C@H](C=O)Cc1c[nH]c2ccccc12",
    "TYR": "N[C@H](C=O)Cc1ccc(O)cc1",
    "HOH": "O",
}


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


def assign_bo_with_template_smiles(
    mol: Chem.Mol, aa_name: str, slice_ids: list, connect_sites: list, max_match_num=10000,
):
    aa_smiles = AA2SMILES[aa_name]
    if isinstance(aa_smiles, str):
        aa_smiles_list = [aa_smiles]
    else:
        aa_smiles_list = aa_smiles

    for idx, aa_smi_ in enumerate(aa_smiles_list):
        params = Chem.SmilesParserParams()
        if "[H]" in aa_smi_:
            params.removeHs = False

        aa_mol = Chem.MolFromSmiles(aa_smi_, params)
        Chem.Kekulize(aa_mol, clearAromaticFlags=True)
        aa_mol2 = Chem.Mol(aa_mol)
        aa_mol2_chg = dict()

        for b in aa_mol2.GetBonds():
            if b.GetBondType() != Chem.BondType.SINGLE:
                b.SetBondType(Chem.BondType.SINGLE)
                b.SetIsAromatic(False)

        # set atom charges to zero;
        for atom in aa_mol2.GetAtoms():
            aa_mol2_chg[atom.GetIdx()] = atom.GetFormalCharge()
            atom.SetFormalCharge(0)

        matches = mol.GetSubstructMatches(aa_mol2, maxMatches=max_match_num)
        filtered_matches = []
        for match_ in matches:
            # if match_[0] in visited_ids: continue
            # visited_ids.update(match_)
            atom = mol.GetAtomWithIdx(match_[0])
            if atom.GetPDBResidueInfo().GetResidueName() == aa_name:
                filtered_matches.append(match_)
    
        if len(filtered_matches) == 0:
            print(f"{aa_smi_}: {idx} not exist!")

        aa_ids = list(range(aa_mol2.GetNumAtoms()))
        for slice_atoms in slice_ids:
            for match_ in filtered_matches:
                if len(set(match_) - set(slice_atoms)) == 0:
                    mapping = dict()
                    for idx1, idx2 in zip(match_, aa_ids):
                        atom = mol.GetAtomWithIdx(idx1)
                        if aa_mol2_chg[idx2] != 0:
                            atom.SetFormalCharge(aa_mol2_chg[idx2])
                        else:
                            if (
                                atom.GetSymbol() == "N"
                                and len(list(atom.GetNeighbors())) > 3
                            ):
                                connect_sites.append(atom.GetIdx())
                        mapping[idx2] = idx1

                    for bond in aa_mol.GetBonds():
                        at1 = bond.GetBeginAtomIdx()
                        at2 = bond.GetEndAtomIdx()
                        new_bond = mol.GetBondBetweenAtoms(mapping[at1], mapping[at2])
                        new_bond.SetBondType(bond.GetBondType())
    return mol


def rdmol_to_omm(rdmol: Chem.Mol) -> app.Modeller:
    # convert RDKit to OpenFF
    from openff.toolkit import Molecule

    off_mol = Molecule.from_rdkit(rdmol, hydrogens_are_explicit=True)

    # convert from OpenFF to OpenMM
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0]
    # convert units from Ångström to nanometers
    mol_positions = mol_positions.to("nanometers")
    # combine topology and positions in modeller object
    omm_mol = app.Modeller(mol_topology, mol_positions)
    return omm_mol


def omm_protein_to_rdmol(
    topology: app.Topology, positions: unit.Quantity = None
) -> Chem.Mol:
    new_mol = Chem.RWMol()
    atom_mapper = dict()
    residue_mapper = defaultdict(list)
    pos_chg_idxs = []
    middle_xt_idxs = []
    for chain in topology.chains():
        chain: app.Chain
        num_res = len(chain._residues)
        for res_idx, residue in enumerate(chain.residues()):
            for atom in residue.atoms():
                atom: app.Atom
                rdatom = Chem.Atom(atom.element.atomic_number)
                # rdatom.SetNoImplicit(True)
                mi = Chem.AtomPDBResidueInfo()
                mi.SetResidueName(atom.residue.name)
                mi.SetResidueNumber(int(atom.residue.id))
                mi.SetChainId(atom.residue.chain.id)
                mi.SetName(format_4letter(atom.name))
                mi.SetInsertionCode("<>")

                if (
                    res_idx == 0 and atom.name == "N"
                ):  # 如果第一个氨基酸不是capped residue, 则设置为N+
                    # 每条链第一个残基的N设置为+1
                    rdatom.SetFormalCharge(1)

                elif (atom.name == "OXT") and (res_idx == num_res - 1):
                    # 末端原子设置为N负
                    rdatom.SetFormalCharge(-1)
                else:
                    rdatom.SetFormalCharge(0)

                rdatom.SetMonomerInfo(mi)
                index = new_mol.AddAtom(rdatom)
                atom_mapper[atom.index] = index
                key = (atom.residue.name, int(atom.residue.id), atom.residue.chain.id)

                if atom.name[-2:] != "XT":  # 封端原子
                    residue_mapper[key].append(index)
                else:
                    if res_idx < num_res - 1:
                        middle_xt_idxs.append(index)

                if res_idx == 1 and atom.name == "N":
                    # 每条链第2个残基的主链N
                    pos_chg_idxs.append(index)

    split_ids = []
    for bond in topology.bonds():
        if bond[0].index in atom_mapper and bond[1].index in atom_mapper:
            at1 = atom_mapper[bond[0].index]
            at2 = atom_mapper[bond[1].index]
            new_mol.AddBond(at1, at2, Chem.BondType.SINGLE)
            if bond[0].residue.id != bond[1].residue.id:
                rdbond = new_mol.GetBondBetweenAtoms(at1, at2)
                split_ids.append(rdbond.GetIdx())

    residue_byres = defaultdict(list)
    for res_info, res_atoms in residue_mapper.items():
        residue_byres[res_info[0]].append(res_atoms)

    connect_sites = []
    # visited_ids = set()
    for res_name, res_atoms_list in residue_byres.items():
        # if res_name=="ACE":
        #     print('ACE')
        assign_bo_with_template_smiles(new_mol, res_name, res_atoms_list, connect_sites)

    # 移除第二个氨基酸的N上多连接的H
    remove_h_idxs = []
    for atidx in list(set(pos_chg_idxs + connect_sites)):
        atom = new_mol.GetAtomWithIdx(atidx)
        count = 0
        tmp_idxs = []
        for nbr_atom in atom.GetNeighbors():
            count += 1
            if nbr_atom.GetSymbol() == "H":
                tmp_idxs.append(nbr_atom.GetIdx())
        remove_h_idxs.extend(tmp_idxs[-(count - 3) :])

    atom_positions = np.array(positions._value) * 10.0
    remove_h_idxs += middle_xt_idxs
    if len(remove_h_idxs) > 0:
        remove_h_idxs = list(sorted(list(set(remove_h_idxs)), reverse=True))
        [new_mol.RemoveAtom(idx) for idx in remove_h_idxs]
        atom_positions = np.delete(atom_positions, remove_h_idxs, axis=0)

    new_mol = new_mol.GetMol()
    new_mol.UpdatePropertyCache(strict=False)
    problems = Chem.DetectChemistryProblems(new_mol)
    for problem in problems:
        cur_idx = problem.GetAtomIdx()
        atom = new_mol.GetAtomWithIdx(cur_idx)
        pdb_info = atom.GetPDBResidueInfo()
        print(
            pdb_info.GetResidueName(),
            pdb_info.GetResidueNumber(),
            pdb_info.GetChainId(),
        )

    Chem.SanitizeMol(new_mol)

    conf = Chem.Conformer(new_mol.GetNumAtoms())
    for aidx in range(new_mol.GetNumAtoms()):
        conf.SetAtomPosition(aidx, atom_positions[aidx])

    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conf, assignId=True)
    new_mol_h = Chem.AddHs(new_mol, addCoords=True)
    # new_mol_h = new_mol

    h_records = defaultdict(int)
    for i in range(new_mol.GetNumAtoms(), new_mol_h.GetNumAtoms()):
        atom = new_mol_h.GetAtomWithIdx(i)
        nbr_atom = atom.GetNeighbors()[0]
        nbr_mi = nbr_atom.GetPDBResidueInfo()
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(nbr_mi.GetResidueName())
        mi.SetResidueNumber(int(nbr_mi.GetResidueNumber()))
        mi.SetChainId(nbr_mi.GetChainId())
        mi.SetInsertionCode("<>")
        nbr_name = nbr_mi.GetName().strip()
        key = (
            nbr_mi.GetResidueName(),
            int(nbr_mi.GetResidueNumber()),
            nbr_mi.GetChainId(),
            nbr_name,
        )
        h_records[key] += 1
        if nbr_name == "N" and nbr_mi.GetResidueName() in ["NME"]:
            label = f"H"
        elif len(nbr_name) == 2:
            label = f"H{nbr_name[1]}{h_records[key]}"
        elif len(nbr_name) == 3:
            label = f"H{nbr_name[1:]}{h_records[key]}"
        else:
            label = f"H{h_records[key]}"
        mi.SetName(format_4letter(label))
        atom.SetMonomerInfo(mi)

    # 重排 H
    residue_mapper = defaultdict(list)
    for atom in new_mol_h.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        key = (
            pdb_info.GetResidueName(),
            pdb_info.GetResidueNumber(),
            pdb_info.GetChainId(),
        )
        residue_mapper[key].append(atom.GetIdx())

    resort_list = []
    [resort_list.extend(values) for _, values in residue_mapper.items()]
    new_mol_h = Chem.RenumberAtoms(new_mol_h, resort_list)
    # Chem.AssignStereochemistryFrom3D(new_mol_h)
    # Chem.SetDoubleBondNeighborDirections(new_mol_h)
    return new_mol_h
