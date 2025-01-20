import os
import argparse
import numpy as np
from dataclasses import asdict
from collections import defaultdict
from posex.utils import DownloadConfig
from posex.data import DatasetGenerator
from posex.preprocess import DataPreprocessor



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, required=True, help="dataset name")
    parser.add_argument("--pdb_ccd_path", type=str, default="/home/lixinze/github/collab/protein_ligand_docking_benchmark/data/dataset/astex_diverse_set_ids.txt")
    parser.add_argument("--download_dir", type=str, default="/data/dataset/posex", help="folder to save the downloaded files")
    args = parser.parse_args()
    np.random.seed(42)
    output_dir = os.path.abspath(args.name)
    assert not os.path.exists(output_dir), f"The {args.name} dataset already exists"
    download_config = asdict(DownloadConfig(args.download_dir))
    pdb_ccd_dict = defaultdict(set)
    with open(args.pdb_ccd_path, "r") as f:
        for pdb_ccd in f.readlines():
            pdb, ccd = pdb_ccd.strip().split("_")
            pdb_ccd_dict[pdb].add(ccd)
    pdbid_list = list(pdb_ccd_dict.keys())
    data_preprocessor = DataPreprocessor(pdbid_list, **download_config)
    pdb_ccd_instance_map = data_preprocessor.run()
    dataset_generator = DatasetGenerator(mode="self_dock",
                                         pdb_ccd_instance_map=pdb_ccd_instance_map,
                                         output_dir=output_dir,
                                         mmseqs_exec=None,
                                         **download_config)
    dataset_generator.select_single_conformation()
    dataset_generator.set_pdb_ccd_dict(pdb_ccd_dict)
    dataset_generator.save_self_dock_res()
