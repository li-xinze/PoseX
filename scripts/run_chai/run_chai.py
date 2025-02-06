import argparse
from pathlib import Path

import torch
from chai_lab.chai1 import run_inference


def main(args: argparse.Namespace):
    run_inference(
        fasta_file=Path(args.fasta_file),
        output_dir=Path(args.output_dir),
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=42,
        device=torch.device(f"cuda:{args.gpu_id}"),
        use_esm_embeddings=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to the fasta file")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to save the output")
    parser.add_argument("--gpu_id", type=int, required=True, help="GPU ID")
    args = parser.parse_args()
    
    main(args)
