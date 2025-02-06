import os
import argparse


def main(args: argparse.Namespace):
    for pdb_ccd in os.listdir(args.input_folder):
        pdb_ccd_dir = os.path.join(args.input_folder, pdb_ccd)
        protein_path = os.path.join(pdb_ccd_dir, f"{pdb_ccd}_protein.pdb")
        ligand_path = os.path.join(pdb_ccd_dir, f"{pdb_ccd}_ligand_start_conf.sdf")
        output_dir = os.path.join(args.output_folder, pdb_ccd)
        os.makedirs(output_dir, exist_ok=True)
        if os.path.isdir(pdb_ccd_dir):
            cmd = f"neuralplexer-inference --task=batched_structure_sampling \
                       --input-receptor {protein_path} \
                       --input-ligand {ligand_path} \
                       --use-template  --input-template {protein_path} \
                       --out-path {output_dir}\
                       --model-checkpoint {args.model_checkpoint} \
                       --n-samples {args.n_samples} \
                       --chunk-size {args.chunk_size} \
                       --num-steps={args.num_steps} \
                       --cuda \
                       --sampler=langevin_simulated_annealing"
            os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_folder", type=str, required=True, help="Path to the input folder")
    parser.add_argument("--output_folder", type=str, required=True, help="Path to the output folder")
    parser.add_argument("--model_checkpoint", type=str, required=True, help="Path to the model checkpoint")
    parser.add_argument("--n_samples", type=int, default=16, help="The number of conformations to generate in total")
    parser.add_argument("--chunk_size", type=int, default=4, help="The number of conformation to generate in parallel")
    parser.add_argument("--num_steps", type=int, default=40,
                        help="The number of steps for the diffusion part of the sampling process")
    args = parser.parse_args()

    main(args)
