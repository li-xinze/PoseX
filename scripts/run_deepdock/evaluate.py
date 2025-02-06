import argparse

from deepdock.prepare_target.computeTargetMesh import compute_inp_surface
from rdkit import Chem
from deepdock.models import *
from deepdock.DockingFunction import dock_compound, get_random_conformation

import numpy as np
import torch


def main(args: argparse.Namespace):
    np.random.seed(123)
    torch.cuda.manual_seed_all(123)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    ligand_model = LigandNet(28, residual_layers=10, dropout_rate=0.10)
    target_model = TargetNet(4, residual_layers=10, dropout_rate=0.10)
    model = DeepDock(ligand_model, target_model, hidden_dim=64, n_gaussians=10, dropout_rate=0.10, dist_threhold=7.).to(
        device)
    checkpoint = torch.load('/DeepDock/Trained_models/DeepDock_pdbbindv2019_13K_minTestLoss.chk',
                            map_location=torch.device(device))
    model.load_state_dict(checkpoint['model_state_dict'])

    ligand_filename = f'{args.pdb_ccd_id}_ligand.mol2'
    sdf_filename = f'{args.pdb_ccd_id}_ligand.sdf'
    target_filename = f'{args.pdb_ccd_id}_protein.pdb'
    target_ply = f'{args.pdb_ccd_id}_protein.ply'
    output_filename = f'{args.pdb_ccd_id}_optimal.sdf'

    compute_inp_surface(target_filename, ligand_filename, dist_threshold=10)

    real_mol = Chem.MolFromMolFile(sdf_filename)
    opt_mol, init_mol, result = dock_compound(real_mol, target_ply, model, dist_threshold=3., popsize=150, seed=123,
                                              device=device)

    writer = Chem.SDWriter(output_filename)
    writer.write(opt_mol, confId=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_ccd_id", type=str, required=True, help="The PDB-CCD id")
    args = parser.parse_args()

    main(args)
