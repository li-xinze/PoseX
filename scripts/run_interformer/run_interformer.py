import os
import os.path
from collections import defaultdict

from rdkit import Chem
import sys
import pandas as pd
import argparse

os.environ['CUDA_VISIBLE_DEVICES'] = "1"


def gen_demo_dock_csv(sdf_f, target, csv_path, isuff=True):
    data = []
    suppl = Chem.SDMolSupplier(sdf_f, sanitize=False)
    for i, mol in enumerate(suppl):
        m_id = mol.GetProp('_Name')
        if mol is not None:
            if isuff:
                data.append([target, 0, i, m_id])
            else:
                data.append([target, i, 0, m_id])

    df = pd.DataFrame(data, columns=['Target', 'pose_rank', 'uff_pose_rank', 'Molecule ID'])
    df.to_csv(csv_path, index=False)


def main(args: argparse.Namespace):
    os.chdir(os.path.join(args.interformer_exec_dir))
    for pdb_ccd in os.listdir(args.input_dir):
        pdb = pdb_ccd.split("_")[0]
        task_dir = os.path.join(args.input_dir, pdb_ccd)
        raw_folder = os.path.join(task_dir, "raw")
        ligand_folder = os.path.join(task_dir, "ligand")
        crystal_ligand_folder = os.path.join(task_dir, "crystal_ligand")
        uff_folder = os.path.join(task_dir, "uff")
        pocket_folder = os.path.join(task_dir, "pocket")
        raw_pocket_folder = os.path.join(raw_folder, "pocket")
        os.makedirs(ligand_folder, exist_ok=True)
        os.makedirs(crystal_ligand_folder, exist_ok=True)
        os.makedirs(uff_folder, exist_ok=True)
        os.makedirs(pocket_folder, exist_ok=True)
        input_protein = os.path.join(raw_folder, f"{pdb_ccd}_protein.pdb")
        input_ligand = os.path.join(raw_folder, f"{pdb_ccd}_ligand.sdf")
        start_conf = os.path.join(raw_folder, f"{pdb_ccd}_ligand_start_conf.sdf")

        crystal_ligand_output_path = os.path.join(crystal_ligand_folder, f"{pdb}_docked.sdf")
        os.system(f"obabel {start_conf} -p 7.4 -O {crystal_ligand_output_path}")
        os.system(f"python tools/rdkit_ETKDG_3d_gen.py {crystal_ligand_folder} {uff_folder}")

        ligand_with_hydrogen_path = os.path.join(ligand_folder, f"{pdb}_docked.sdf")
        os.system(f"obabel {input_ligand} -p 7.4 -O {ligand_with_hydrogen_path}")
        # os.system(f"python tools/rdkit_ETKDG_3d_gen.py {ligand_folder} {uff_folder}")

        reduced_protein_path = os.path.join(raw_pocket_folder, f"{pdb}_reduce.pdb")
        os.system(f"mkdir -p {raw_pocket_folder} && reduce {input_protein} > {reduced_protein_path}")

        pocket_output_path = os.path.join(raw_pocket_folder, "output", f"{pdb}_pocket.pdb")
        os.system(
            f"python tools/extract_pocket_by_ligand.py {raw_pocket_folder} {ligand_folder} 0 "
            f"&& mv {pocket_output_path} {pocket_folder}")

        energy_output = os.path.join(task_dir, "energy_output")
        csv_path = os.path.join(task_dir, "demo_dock.csv")

        gen_demo_dock_csv(input_ligand, pdb, csv_path)

        predict_energy_cmd = f"PYTHONPATH=interformer/ python inference.py -test_csv {csv_path} \
                            -work_path {task_dir} \
                            -ensemble checkpoints/v0.2_energy_model \
                            -batch_size 1 \
                            -posfix *val_loss* \
                            -energy_output_folder {energy_output} \
                            -reload \
                            -debug"
        os.system(predict_energy_cmd)
        # os.system(f"cp {start_conf} {energy_output}/uff/{pdb}_uff.sdf")
        os.system(
            f'OMP_NUM_THREADS="64,64" python docking/reconstruct_ligands.py -y --cwd {energy_output} -y --find_all --uff_folder uff find')
        os.system(f'python docking/reconstruct_ligands.py --cwd {energy_output} --find_all stat')
        os.system(
            f'python docking/merge_summary_input.py {os.path.join(energy_output, "ligand_reconstructing/stat_concated.csv")} {csv_path}')
        infer_dir = os.path.join(task_dir, "infer")
        os.makedirs(infer_dir, exist_ok=True)
        os.system(f'cp -r {os.path.join(energy_output, "ligand_reconstructing")} {infer_dir}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True, help="Path to the input files")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--interformer_exec_dir", type=str, required=True, help="Path to the Unimol codebase")
    args = parser.parse_args()

    main(args)
