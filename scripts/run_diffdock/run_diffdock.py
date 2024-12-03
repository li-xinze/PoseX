import argparse
import subprocess


def main(args: argparse.Namespace):

    subprocess.run(
        [
            "python",
            "inference.py",
            "--protein_ligand_csv",
            args.input_csv_path,
            "--out_dir",
            args.output_dir,
            "--cuda_device_index",
            str(args.gpu_id),
            "--model_dir",
            args.model_dir,
            "--confidence_model_dir",
            args.confidence_model_dir, 
            "--inference_steps",
            "20",
            "--samples_per_complex",
            "40",
            "--actual_steps",
            "18",
            "--no_final_step_noise",
        ],
        cwd=args.diffdock_exec_dir,
        check=True
    )  # nosec



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv_path", type=str, required=True, help="Path to the protein_ligand_csv file")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--diffdock_exec_dir", type=str, required=True, help="Path to the DiffDock codebase")
    parser.add_argument("--model_dir", type=str, required=True, help="Path to the model_dir")
    parser.add_argument("--confidence_model_dir", type=str, required=True, help="Path to the confidence_model_dir")
    parser.add_argument("--gpu_id", type=int, required=True, help="GPU ID")
    args = parser.parse_args()
    
    main(args)