import os
import argparse
import numpy as np
from dataclasses import asdict
from posex.utils import DownloadConfig
from posex.data import DatasetGenerator
from posex.preprocess import DataPreprocessor

        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str, default="cross_dock", help="benchmark mode(self_dock or cross_dock)")
    parser.add_argument("--pdbid_path", type=str, default="/data/dataset/posex/pdbid_cross_2.0.txt")
    parser.add_argument("--download_dir", type=str, default="/data/dataset/posex", help="folder to save the downloaded files")
    parser.add_argument("--mmseqs_exec", type=str, default="/data/software/mmseqs2/mmseqs/bin/mmseqs", help="path to mmseqs exec")
    args = parser.parse_args()
    np.random.seed(42)
    output_dir = os.path.abspath(args.mode)
    assert not os.path.exists(output_dir), f"The {args.mode} dataset already exists"
    download_config = asdict(DownloadConfig(args.download_dir))
    with open(args.pdbid_path, "r") as f:
        pdbid_list = f.readline().strip().split(",")
    data_preprocessor = DataPreprocessor(pdbid_list, **download_config)
    pdb_ccd_instance_map = data_preprocessor.run()
    dataset_generator = DatasetGenerator(mode=args.mode,
                                         pdb_ccd_instance_map=pdb_ccd_instance_map,
                                         output_dir=output_dir,
                                         mmseqs_exec=args.mmseqs_exec,
                                         **download_config)
    dataset_generator.run()
