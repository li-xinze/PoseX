import os
import argparse


def main(args: argparse.Namespace):
    os.chdir(os.path.join(args.interformer_exec_dir))
    for pdb_ccd in os.listdir(args.input_dir):
        task_dir = os.path.join(args.input_dir, pdb_ccd)
        raw_folder = os.path.join(task_dir, "raw")
        ligand_folder = os.path.join(task_dir, "ligand")
        utf_folder = os.path.join(task_dir, "utf")
        pocket_folder = os.path.join(task_dir, "pocket")
        raw_pocket_folder = os.path.join(raw_folder, "pocket")
        os.makedirs(ligand_folder, exist_ok=True)
        os.makedirs(utf_folder, exist_ok=True)
        os.makedirs(pocket_folder, exist_ok=True)
        input_protein = os.path.join(raw_folder, f"{pdb_ccd}_protein.pdb")
        input_ligand = os.path.join(raw_folder, f"{pdb_ccd}_ligand/_start_conf.sdf")

        ligand_with_hydrogen_path = os.path.join(ligand_folder, f"{pdb_ccd}_docked.sdf")
        os.system(f"obabel {input_ligand} -p 7.4 -O {ligand_with_hydrogen_path}")
        os.system(f"python tools/rdkit_ETKDG_3d_gen.py {ligand_folder} {utf_folder}")

        reduced_protein_path = os.path.join(raw_pocket_folder, f"{pdb_ccd}_reduce.pdb")
        os.system(f"mkdir -p {raw_pocket_folder} && reduce {input_protein} > {reduced_protein_path}")

        pocket_output_path = os.path.join(raw_pocket_folder, "output", f"{pdb_ccd}_pocket.pdb")
        os.system(
            f"python tools/extract_pocket_by_ligand.py {raw_pocket_folder} {ligand_folder} 1 "
            f"&& mv {pocket_output_path} {pocket_folder}")

        energy_output = os.path.join(task_dir, "energy_output")
        csv_path = os.path.join(task_dir, "demo_dock.csv")
        predict_energy_cmd = f"python inference.py -test_csv {csv_path} \
                            -work_path examples/ \
                            -ensemble checkpoints/v0.2_energy_model \
                            -batch_size 1 \
                            -posfix *val_loss* \
                            -energy_output_folder {energy_output} \
                            -reload \
                            -debug"
        os.system(predict_energy_cmd)

        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True, help="Path to the input files")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--interformer_exec_dir", type=str, required=True, help="Path to the Unimol codebase")
    parser.add_argument("--ckpt_path", type=str, required=True, help="Path to the model chekpoint")
    args = parser.parse_args()

    main(args)
