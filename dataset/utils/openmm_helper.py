import argparse
import openmm as omm
import openmm.app as app

from pathlib import Path
import pdbfixer
from rdkit import Chem
from openff.toolkit.topology import Molecule
from collections import defaultdict
import numpy as np
import copy,re
from openmm import unit, Vec3
from openmm.app import PDBFile
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
    SMIRNOFFTemplateGenerator,
)
from io import StringIO
from functools import partial
from subprocess import check_output as _call

from typing import *

from dataset.utils.common_helper import create_logger, temprary_filename

logger = create_logger(__name__)

METALS = ["Na", "K", "Ca", "Mg", "Fe", "Zn", "Cu", "Mn", "Co", "Ni"]
run_shell_cmd = partial(_call, shell=True)


def assign_mol_with_pos(mol: Chem.Mol, pos: np.array, is_pdb=False) -> Chem.Mol:
    conf = mol.GetConformer()
    conf.SetPositions(pos)
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conf, assignId=True)
    return new_mol

def load_sdf_to_omm(
    sdf_fn: str | Path,
) -> app.Modeller:
    # convert RDKit to OpenFF
    rdk_mol = Chem.SDMolSupplier(sdf_fn, removeHs=False)[0]
    rdk_mol = Chem.AddHs(rdk_mol, addCoords=True)
    off_mol = Molecule.from_rdkit(
        rdk_mol, hydrogens_are_explicit=True, allow_undefined_stereo=True
    )
    props = rdk_mol.GetPropsAsDict()
    [rdk_mol.ClearProp(k) for k in props.keys()]

    # convert from OpenFF to OpenMM
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0]
    # convert units from Ångström to nanometers
    mol_positions = mol_positions.to("nanometers")
    lig_pos = [Vec3(*x.m) for x in mol_positions] * unit.nanometers
    # combine topology and positions in modeller object
    for residue in mol_topology.residues():
        residue.name = "UNL"
    omm_mol = app.Modeller(mol_topology, lig_pos)
    off_mol.name = 'UNL'
    return omm_mol, off_mol, rdk_mol


def get_am1bcc_charge(ligand_fn: Path):
    lig_mol = Chem.SDMolSupplier(ligand_fn, removeHs=False)[0]
    total_charge = 0
    for atom in lig_mol.GetAtoms():
        total_charge += atom.GetFormalCharge()
    with temprary_filename(mode='w', suffix='_ligand.mol2') as tmp_out_fn:
        cmd1 =f'antechamber -i {ligand_fn} -fi sdf -o {tmp_out_fn} -fo mol2 -c bcc -nc {total_charge} -at gaff2 -ek "qm_theory=\'AM1\', scfconv=1.d-10, maxcyc=0, grms_tol=0.0005, ndiis_attempts=700" -pf y'
        run_shell_cmd(cmd1)
        with open(tmp_out_fn, 'r') as f:
            mol_data = f.read()
        atom_section = re.search(r'@<TRIPOS>ATOM([\s\S]+?)@<TRIPOS>BOND', mol_data)
        atom_lines = atom_section.group(1).strip().split("\n")
        partial_charges = [float(line.split()[-1]) for line in atom_lines]

    assert len(partial_charges) == lig_mol.GetNumAtoms(), f"{ligand_fn} partial charges are not correct"
    
    from openff.toolkit import Quantity
    from openff.toolkit import  unit as off_unit
    partial_charges = Quantity(np.asarray(partial_charges), off_unit.elementary_charge)
    return partial_charges



class ProLigRelax:
    def __init__(
        self,
        protein_mol: Chem.Mol = None,
        platform="CPU:16",
        ligand_ff="gaff",
        receptor_ff="amber",
        charge_name="mmff94",
        missing_residues: list = None,
        is_constrain: bool = False
    ) -> None:
        self.receptor_ffname = receptor_ff
        self.ligand_ffname = ligand_ff
        self.partial_chargename = charge_name

        self.platform_properties = {}
        if platform.startswith("CPU"):
            device, num_core = platform.split(":")
            self.platform = omm.Platform.getPlatformByName(device)
            self.platform.setPropertyDefaultValue("Threads", str(num_core))

        elif platform.startswith("CUDA"):
            device, device_id = platform.split(":")
            self.platform = omm.Platform.getPlatformByName(device)
            self.platform_properties.update({"DeviceIndex": device_id})
            self.platform_properties.update({"Precision": "single"})

        else:
            raise NotImplementedError()

        self.receptor_rdmol = protein_mol
        self.receptor_omm = app.PDBFile(
            StringIO(Chem.MolToPDBBlock(protein_mol).replace("< ", "  "))
        )
        self.base_forcefield = self._load_rec_forcefield()
        self.missing_residues = missing_residues
        self.is_constrain = is_constrain
        self.forcefield_kwargs = {
            "nonbondedMethod": app.NoCutoff,
            "nonbondedCutoff": 2.0 * unit.nanometer,  # (default: 1)
            "rigidWater": False,
            "removeCMMotion": True,
            "constraints": app.HBonds,
            "soluteDielectric": 1.0,
            "solventDielectric": 80.0,
        }

    def _load_rec_forcefield(self) -> None:
        if self.receptor_ffname == "amber":
            receptor_ffs = [
                "amber14/protein.ff14SB.xml", 
                "amber14/tip3pfb.xml", 
                "implicit/obc2.xml"
            ]
            forcefield = app.ForceField(*receptor_ffs)
        else:
            raise NotImplementedError()

        return forcefield

    def _add_lig_forcefield(self, ligand_mol: Molecule, ligand_fn: Path = None) -> app.ForceField:
        if getattr(ligand_mol, "partial_charges") is None:
            if self.partial_chargename != "am1bcc":
                ligand_mol.assign_partial_charges(
                    partial_charge_method=self.partial_chargename,
                    use_conformers=ligand_mol.conformers[0],
                )
            else:
                # ligand_mol.assign_partial_charges(partial_charge_method='mmff94')
                ligand_mol.partial_charges=get_am1bcc_charge(ligand_fn)

        if self.ligand_ffname == "gaff":
            ffgen = GAFFTemplateGenerator(forcefield="gaff-2.11")
            ffxml_contents = ffgen.generate_residue_template(ligand_mol)

        elif self.ligand_ffname == "openff":
            ffgen = SMIRNOFFTemplateGenerator(forcefield="openff-2.1.0")
            ffxml_contents = ffgen.generate_residue_template(ligand_mol)

        else:
            raise NotImplementedError()

        forcefield: app.ForceField = copy.deepcopy(self.base_forcefield)
        forcefield.loadFile(StringIO(ffxml_contents))
        return forcefield

    @staticmethod
    def constrain_pocket(
        system, top, missing_residues: list = [], level="main"
    ) -> None:
        assert level in ["none", "main", "heavy", "all"]
        for atom in top.topology.atoms():
            if atom.residue.name in ["UNL", "HOH"] + METALS:  # water和metal不被限制
                continue

            if len(missing_residues) > 0:  # missing residues不设置限制
                tag = (atom.residue.chain.index, atom.residue.id, atom.residue.name)
                if tag in missing_residues:
                    continue

            if level == "all":
                system.setParticleMass(atom.index, 0)
            elif level == "main" and atom.name in ["CA", "C", "N"]:
                system.setParticleMass(atom.index, 0)
            elif level == "heavy" and atom.name[0] != "H":
                system.setParticleMass(atom.index, 0)

    @staticmethod
    def restrain_pocket(
        system,
        model: app.Modeller,
        missing_residues: list = [],
        level="main",
        stiffness: float = 10,
    ) -> None:
        assert level in ["none", "main", "heavy", "all"]

        forces = omm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
        forces.addGlobalParameter("k", stiffness)
        for p in ["x0", "y0", "z0"]:
            forces.addPerParticleParameter(p)

        for i, atom in enumerate(model.topology.atoms()):
            if atom.residue.name in ["UNL", "HOH"] + METALS:
                continue

            if len(missing_residues) > 0:
                tag = (atom.residue.chain.index, atom.residue.id, atom.residue.name)
                if tag in missing_residues:
                    continue

            if level == "all":
                forces.addParticle(i, model.positions[i])
            elif level == "main" and atom.name in ["CA", "C", "N"]:
                forces.addParticle(i, model.positions[i])
            elif level == "heavy" and atom.name[0] != "H":
                forces.addParticle(i, model.positions[i])

        system.addForce(forces)

    def _minimize_energy(
        self, cplx_forcefield: app.ForceField, cplx_model: app.Modeller) -> np.ndarray:    
        if self.is_constrain:
            self.forcefield_kwargs["constraints"]= None
            system = cplx_forcefield.createSystem(
                cplx_model.topology, **self.forcefield_kwargs
            )
            self.constrain_pocket(
                system, cplx_model, missing_residues=self.missing_residues, 
                level="main"
            )
        else:
            system = cplx_forcefield.createSystem(
                cplx_model.topology, **self.forcefield_kwargs
            )
            self.restrain_pocket(
                system, cplx_model, missing_residues=self.missing_residues, 
                level="main"
            )


        integrator = omm.LangevinIntegrator(
            300, 1, 0.001
        )  # from heyi mononor parameters
        # only use one cpu per relaxation
        simulation = app.Simulation(
            cplx_model.topology,
            system,
            integrator,
            self.platform,
            self.platform_properties,
        )
        simulation.context.setPositions(cplx_model.positions)
        logger.info(
            f"Minimizing with {self.receptor_ffname} (protein) + {self.ligand_ffname} (ligand) + {self.partial_chargename} charge..."
        )
        simulation.minimizeEnergy()

        state = simulation.context.getState(getPositions=True)
        minimized_positions = state.getPositions(asNumpy=True)
        return minimized_positions

    def prepare_one_cplx_and_relax(self, ligand_fn: Path) -> None:
        ligand_omm, ligand_off, ligand_rdk = load_sdf_to_omm(ligand_fn)
        receptor_rdmol = Chem.Mol(self.receptor_rdmol)

        cplx_model = app.Modeller(
            self.receptor_omm.topology, self.receptor_omm.positions
        )
        cplx_model.add(ligand_omm.topology, ligand_omm.positions)
        logger.info(f"{len(cplx_model.positions)=}")

        cplx_forcefield = self._add_lig_forcefield(ligand_off, ligand_fn)
        relaxed_pos = self._minimize_energy(cplx_forcefield, cplx_model)

        angstroms_pos = relaxed_pos.value_in_unit(unit.angstroms)
        nanometers_pos = relaxed_pos.value_in_unit(unit.nanometers)

        rec_num = len(self.receptor_omm.positions)
        pos1 = angstroms_pos[:rec_num]
        pos2 = angstroms_pos[rec_num:]

        self.receptor_omm.positions = nanometers_pos[:rec_num] * unit.nanometers

        protein_mol = assign_mol_with_pos(receptor_rdmol, pos1, is_pdb=True)
        ligand_mol = assign_mol_with_pos(ligand_rdk, pos2)

        return protein_mol, ligand_mol

    def run_batch(self, batch_fns: List[Path]) -> None:
        for ligand_ in batch_fns:
            self.prepare_one_cplx_and_relax(ligand_)

        return f"finished: {batch_fns}"
