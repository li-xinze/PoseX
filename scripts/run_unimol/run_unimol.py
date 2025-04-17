import os
import argparse



def main(args: argparse.Namespace):
    os.chdir(os.path.join(args.unimol_exec_dir, "interface"))
    for pdb_ccd in os.listdir(args.input_dir):
        input_protein = os.path.join(args.input_dir, pdb_ccd, f"{pdb_ccd}_protein.pdb")
        input_docking_grid = os.path.join(args.input_dir, pdb_ccd, f"{pdb_ccd}.json")
        input_ligand = os.path.join(args.input_dir, pdb_ccd, f"{pdb_ccd}_ligand_start_conf.sdf")
        output_dir = os.path.join(args.output_dir, pdb_ccd)
        os.makedirs(output_dir, exist_ok=True)
        cmd = f"python demo.py --mode single --conf-size 10 --cluster \
        --input-protein {input_protein} \
        --input-ligand {input_ligand} \
        --input-docking-grid {input_docking_grid} \
        --output-ligand-name ligand_predict \
        --output-ligand-dir {output_dir} \
        --steric-clash-fix \
        --model-dir {args.ckpt_path}"
        os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True, help="Path to the input files")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--unimol_exec_dir", type=str, required=True, help="Path to the Unimol codebase")
    parser.add_argument("--ckpt_path", type=str, required=True, help="Path to the model chekpoint")
    args = parser.parse_args()

    main(args)
