import os
import argparse


def main(args: argparse.Namespace):
    for pdb_ccd in os.listdir(args.input_folder):
        os.chdir(os.path.join(args.input_folder, pdb_ccd))
        os.system(f"python evaluate.py --pdb_ccd {pdb_ccd}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder")
    args = parser.parse_args()

    main(args)
