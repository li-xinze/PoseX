import os
import argparse


def main(args: argparse.Namespace):
    with open(args.evaluation_template_path) as f:
        data = f.read()

    for pdb_ccd in os.listdir(args.input_folder):
        if os.path.isdir(pdb_ccd):
            new_data = data.replace("LIGAND", pdb_ccd)
            evaluation_output = os.path.join(args.input_dir_path, pdb_ccd, "evaluation.py")
            with open(evaluation_output, "w") as fw:
                fw.write(new_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder")
    parser.add_argument("--evaluation_template_path", type=str, required=True, help="Path to evaluation template")
    args = parser.parse_args()

    main(args)
