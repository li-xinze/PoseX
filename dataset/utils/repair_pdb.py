import argparse
from tracemalloc import start
import openmm as omm
import openmm.app as app

from pathlib import Path
import pdbfixer
from rdkit import Chem
from openff.toolkit.topology import Molecule
from collections import defaultdict
import numpy as np
import copy
from openmm import unit, Vec3
from openmm.app import PDBFile
from openmmforcefields.generators import GAFFTemplateGenerator
from io import StringIO
from typing import *

from dataset.utils.mol_correct_helper import omm_protein_to_rdmol
from dataset.utils.common_helper import create_logger

logger = create_logger(__name__)


class PDBFixer(pdbfixer.PDBFixer):
    def addMissingHydrogens(self, pH=7.0, forcefield=None) -> None:
        extraDefinitions = self._downloadNonstandardDefinitions()
        variants = [
            (
                self._describeVariant(res, extraDefinitions)
                if res.name not in ["ACE", "NME"]
                else None
            )
            for res in self.topology.residues()
        ]
        modeller = app.Modeller(self.topology, self.positions)
        modeller.addHydrogens(pH=pH, forcefield=forcefield, variants=variants)
        self.topology = modeller.topology
        self.positions = modeller.positions


def fix_pdb(pdb_fn: Path, cap_n_ter="ACE", cap_c_ter="NME") -> Chem.Mol:
    fixer = PDBFixer(filename=str(pdb_fn))

    fixer.findNonstandardResidues()
    logger.info(f"{fixer.nonstandardResidues=}")
    fixer.replaceNonstandardResidues()

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    logger.info(f"{fixer.missingResidues=}")
    logger.info(f"{fixer.missingAtoms=}")
    logger.info(f"{fixer.missingTerminals=}")

    for i in range(len(fixer.topology._chains)):
        chain = fixer.topology._chains[i]
        if chain._residues[0].name == "HOH":
            continue

        first_head = (i, 0)
        if first_head not in fixer.missingResidues:
            fixer.missingResidues[first_head] = [cap_n_ter]
        else:
            fixer.missingResidues[first_head] = [cap_n_ter] + fixer.missingResidues[
                first_head
            ][-1:]

        last_head = (i, len(chain))
        if last_head not in fixer.missingResidues:
            fixer.missingResidues[last_head] = [cap_c_ter]
        else:
            fixer.missingResidues[last_head] = fixer.missingResidues[last_head][-1:] + [
                cap_c_ter
            ]

    fixer.missingTerminals = {}
    logger.info(f"{fixer.missingResidues=}")

    start_num = 0
    missing_residues = []
    for chain in fixer.topology.chains():
        if chain._residues[0].name == "HOH":
            continue
        for residue in chain.residues():
            tag = (residue.chain.index, residue.index - start_num)
            if tag in fixer.missingResidues:
                cur_residues = fixer.missingResidues[tag]
                cur_len = len(cur_residues)
                [
                    missing_residues.append(
                        (residue.chain.index, str(int(residue.id) - (cur_len - i)), res)
                    )
                    for i, res in enumerate(cur_residues)
                ]

            if (residue.index + 1 - start_num) == len(chain):
                tag = (residue.chain.index, residue.index + 1 - start_num)
                cur_residues = fixer.missingResidues[tag]
                [
                    missing_residues.append(
                        (residue.chain.index, str(int(residue.id) + (i + 1)), res)
                    )
                    for i, res in enumerate(cur_residues)
                ]
                start_num += len(chain)

    logger.info(f"{len(list(fixer.topology.atoms()))=}")
    fixer.addMissingAtoms(seed=0)
    logger.info(f"add {missing_residues=}")
    logger.info(f"{len(list(fixer.topology.atoms()))=}")

    protein_mol = omm_protein_to_rdmol(fixer.topology, fixer.positions)
    return protein_mol, missing_residues


def test_fix_cif(pdbx_fn: Path, cap_n_ter="ACE", cap_c_ter="NME"):
    with open(pdbx_fn) as f:
        fixer = pdbfixer.PDBFixer(pdbxfile=f)
    # fixer = pdbfixer.PDBFixer(pdbxfile=str(pdbx_fn))
    fixer.findNonstandardResidues()
    logger.info(f"{fixer.nonstandardResidues=}")
    # fixer.replaceNonstandardResidues()

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    logger.info(f"{fixer.missingResidues=}")
    logger.info(f"{fixer.missingAtoms=}")
    logger.info(f"{fixer.missingTerminals=}")

    for i in range(len(fixer.topology._chains)):
        chain = fixer.topology._chains[i]
        if chain._residues[0].name == "HOH":
            continue

        first_head = (i, 0)
        if first_head not in fixer.missingResidues:
            fixer.missingResidues[first_head] = [cap_n_ter]
        else:
            fixer.missingResidues[first_head] = [cap_n_ter] + fixer.missingResidues[
                first_head
            ][-1:]

        last_head = (i, len(chain))
        if last_head not in fixer.missingResidues:
            fixer.missingResidues[last_head] = [cap_c_ter]
        else:
            fixer.missingResidues[last_head] = fixer.missingResidues[last_head][-1:] + [
                cap_c_ter
            ]

    fixer.missingTerminals = {}
    logger.info(f"{fixer.missingResidues=}")

    logger.info(f"{len(list(fixer.topology.atoms()))=}")
    fixer.addMissingAtoms(seed=0)
    # logger.info(f'{missing_res}')
    logger.info(f"{len(list(fixer.topology.atoms()))=}")
    # fixer.addMissingHydrogens(7.0)

    protein_mol = omm_protein_to_rdmol(fixer.topology, fixer.positions)

    return protein_mol


if __name__ == "__main__":
    fix_pdb()
    test_fix_cif(
        "/hpc-cache-pfs/home/xyj/code/protein_ligand_docking_benchmark/posex/mmcif_raw/8UCB_X1T.cif"
    )
