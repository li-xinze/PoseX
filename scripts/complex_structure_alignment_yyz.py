from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def sdf_to_pdb(sdf_file, output_dir):
    # Read the SDF file
    suppl = Chem.SDMolSupplier(sdf_file)
    for i,mol in enumerate(suppl):
        pdb_file = os.path.join(output_dir, f'gt_ligand{i+1}.pdb')
        if mol is None:
            raise ValueError("Invalid or empty SDF file")
        
        # Generate 3D coordinates if missing
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol)
        
        # Write the molecule to a PDB file
        with open(pdb_file, "w") as f:
            f.write(Chem.MolToPDBBlock(mol))
    return len(suppl)
def calculate_ligand_rmsd(aligned_native_pdb, aligned_pred_pdb):
    atom_coords_native = {}
    atom_coords_pred = {}
    with open(aligned_native_pdb, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("HETATM"):
                atom_id = line.split()[2]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_coords_native[atom_id] = [x, y, z]
    with open(aligned_pred_pdb, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("HETATM"):
                atom_id = line.split()[2].replace("_", "")
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_coords_pred[atom_id] = [x, y, z]
    coords_native = []
    coords_pred = []
    for atom_key in atom_coords_native:
        if atom_key in atom_coords_pred:
            coords_native.append(atom_coords_native[atom_key])
            coords_pred.append(atom_coords_pred[atom_key])
    rmsd = np.sqrt(np.mean(np.sum((np.array(coords_native) - np.array(coords_pred))**2, axis=1)))
    return rmsd


def align_with_pkt(native_protein, native_ligand, predicted_protein, predicted_ligand, output_dir, radius=10.0, ligand_resname="UNL"):


    # Load the structures
    pkt_id = os.path.split(native_ligand)[1][:-4].replace("gt_ligand", "")
    pymol.finish_launching()
    cmd.load(native_protein, "native_protein")
    cmd.load(predicted_protein, "pred_protein")

    cmd.load(native_ligand_pdb, "native_ligand")
    cmd.load(predicted_ligand_pdb, "predicted_ligand")
    cmd.create("native", "native_protein or native_ligand")  
    cmd.create("predicted", "pred_protein or predicted_ligand")

    cmd.select("gt_ligand", f"native and resn {ligand_resname}")
    cmd.select("pred_ligand", f"predicted and resn LIG")

    
    # Select nearby N, CA, and C atoms within the radius
    cmd.select("native_pkt_atoms", f"(native and name CA+N+C) within {radius} of gt_ligand")
    cmd.select("predicted_pkt_atoms", f"(predicted and name CA+N+C) within {radius} of pred_ligand")

    # Align the nearby atoms
    cmd.align("predicted_pkt_atoms", "native_pkt_atoms")

    aligned_pred_ligand_pdb = os.path.join(output_dir, f'aligned_pred_ligand_{pkt_id}.pdb')
    aligned_gt_ligand_pdb = os.path.join(output_dir, f'aligned_gt_ligand_{pkt_id}.pdb')
    aligned_pred_protein_pdb = os.path.join(output_dir, f'aligned_pred_protein_{pkt_id}.pdb')
    aligned_gt_protein_pdb = os.path.join(output_dir, f'aligned_gt_protein_{pkt_id}.pdb')

    cmd.save(aligned_pred_ligand_pdb, "pred_ligand")
    cmd.save(aligned_gt_ligand_pdb, "gt_ligand")
    cmd.save(aligned_pred_protein_pdb, "predicted")
    cmd.save(aligned_gt_protein_pdb, "native")
    
    # Calculate ligand RMSD
    
    cmd.alter("all", "segi=''")
    cmd.alter("all", "chain=''")
    ligand_rmsd = cmd.rms_cur('pred_ligand', 'gt_ligand', matchmaker=4)
    
    # ligand_rmsd = cmd.rms_cur('pred_ligand', 'gt_ligand', matchmaker=4)
    
    #ligand_rmsd = cmd.align("pred_ligand", "gt_ligand")[0]
    
    cmd.delete("all")
    
    # Clean up PyMOL session
    
    return ligand_rmsd, aligned_pred_ligand_pdb, aligned_gt_ligand_pdb

# Example usage

output_root = 'protein_ligand_docking/output' 
native_root = 'protein_ligand_docking/astex_diverse_set'
chai_pred_root = 'protein_ligand_docking/data/benchmark/astex/chai/output'

targets = os.listdir(chai_pred_root)
rmsd_dict = {}

for target in targets:
    best_rmsd = 1000000
    output_dir = os.path.join(output_root, target)
    os.makedirs(output_dir, exist_ok=True)
    predicted_protein_pdb = os.path.join(chai_pred_root, target, f'{target}_model_protein.pdb')
    predicted_ligand_pdb = os.path.join(chai_pred_root, target, f'{target}_model_ligand.pdb')
    native_protein_pdb = os.path.join(native_root, target, f'{target}_protein.pdb')
    native_ligands_sdf = os.path.join(native_root, target, f'{target}_ligands.sdf')

    mumber_of_ligands = sdf_to_pdb(native_ligands_sdf, output_dir)
    for i in range(mumber_of_ligands):
        native_ligand_pdb = os.path.join(output_dir, f'gt_ligand{i+1}.pdb')
        ligand_rmsd_cur, aligned_pred_ligand_pdb, aligned_gt_ligand_pdb = align_with_pkt(native_protein_pdb, native_ligand_pdb, \
                                                        predicted_protein_pdb, predicted_ligand_pdb, output_dir)
        ligand_rmsd = calculate_ligand_rmsd(aligned_gt_ligand_pdb, aligned_pred_ligand_pdb)
        if ligand_rmsd < best_rmsd:
            best_rmsd = ligand_rmsd
    print(f'{target}: {best_rmsd}')
    rmsd_dict[target] = best_rmsd
with open(os.path.join(output_root,'chai_rmsd.txt'), 'w') as f:
    for target, rmsd in rmsd_dict.items():
        f.write(f'{target}: {rmsd}\n')