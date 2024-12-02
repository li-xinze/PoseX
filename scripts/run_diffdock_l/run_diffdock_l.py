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
            "--config",
            args.config_path,
        ],
        cwd=args.diffdock_exec_dir,
        check=True
    )



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv_path", type=str, required=True, help="Path to the protein_ligand_csv file")
    parser.add_argument("--diffdock_exec_dir", type=str, required=True, help="Path to the DiffDock_L codebase")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--config_path", type=str, required=True, help="Path to the config file")
    args = parser.parse_args()
    
    main(args)