import os
import shutil
import argparse
import subprocess


def main(args: argparse.Namespace):
    os.environ["MKL_THREADING_LAYER"] = "GNU"  
    header = args.itemname
    result_root = os.path.join(args.dynamicbind_exec_dir, "inference", "outputs", "results")
    result_dir = os.path.join(result_root, header)
    subprocess.run(
        [
            "python",
            "run_single_protein_inference.py",
            args.protein_filepath,
            args.ligand_filepath,
            "--samples_per_complex",
            "40",
            "--header",
            header,
            "--device",
            str(args.gpu_id),
            "--python",
            "python",
            "--relax_python",
            "python",
            "--results",
            result_root,
            "--no_relax",
            "--paper", 
        ],
        cwd=args.dynamicbind_exec_dir,
        check=True
    )  

    shutil.move(result_dir, args.output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--itemname", type=str, required=True, help="PDB_CCD_ID")
    parser.add_argument("--protein_filepath", type=str, required=True, help="Path to the protein pdb file")
    parser.add_argument("--ligand_filepath", type=str, required=True, help="Path to the ligand csv file")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--dynamicbind_exec_dir", type=str, required=True, help="Path to the DynamicBind project")
    parser.add_argument("--gpu_id", type=int, required=True, help="GPU ID")
    args = parser.parse_args()
    
    main(args)