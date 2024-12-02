import os
import argparse
import subprocess


def main(args: argparse.Namespace):
    save_pt_dir = os.path.join(args.output_dir, "temp_files")
    save_mols_dir = os.path.join(save_pt_dir, "mol")
    # preprocess ligand
    subprocess.run(
        [
            "python",
            "inference_preprocess_mol_confs.py",
            "--index_csv",
            args.input_csv_path,
            "--save_mols_dir",
            save_mols_dir,
            "--num_threads",
            "1",
        ],
        cwd=args.fabind_exec_dir,
        check=True
    )
    # preprocess protein
    subprocess.run(
        [
            "python",
            "inference_preprocess_protein.py",
            "--pdb_file_dir",
            args.input_data_dir,
            "--save_pt_dir",
            save_pt_dir,
            "--cuda_device_index",
            str(args.gpu_id),
        ],
        cwd=args.fabind_exec_dir,
        check=True
    ) 
    # inference
    subprocess.run(
        [
            "python",
            "fabind_inference.py",
            "--ckpt",
            args.ckpt_path,
            "--batch_size",
            "4",
            "--seed",
            "42",
            "--test-gumbel-soft",
            "--redocking",
            "--post-optim",
            "--write-mol-to-file",
            "--sdf-output-path-post-optim",
            args.output_dir,
            "--index-csv",
            args.input_csv_path,
            "--preprocess-dir",
            save_pt_dir,
            "--cuda_device_index",
            str(args.gpu_id),
        ],
        cwd=args.fabind_exec_dir,
        check=True
    ) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv_path", type=str, required=True, help="Path to the ligand_csv file")
    parser.add_argument("--input_data_dir", type=str, required=True, help="Path to the protein pdb dir")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--fabind_exec_dir", type=str, required=True, help="Path to the FABind codebase")
    parser.add_argument("--ckpt_path", type=str, required=True, help="Path to the model chekpoint")
    parser.add_argument("--gpu_id", type=int, required=True, help="GPU ID")
    args = parser.parse_args()
    
    main(args)
