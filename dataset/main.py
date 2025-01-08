import os
import argparse
import numpy as np
from posex.data import DatasetGenerator
from posex.preprocess import DataPreprocessor
from dataclasses import dataclass, fields, asdict


@dataclass
class DownloadConfig:
    """Download directory configuration

    Attributes
        download_dir: str, folder to save the downloaded files
        ccd_dir : str, default="ccd", subfolder to save ccd files
        ccd_path: str, default="ccd.csv", path to save ccd table
        bcif_dir: str, default="bcif", subfolder to save pdb entries (bcif format)
        cif_dir: str, default="cif", subfolder to save pdb entries (cif format)
        vs_dir: str, default="entry", subfolder to save json files containing validation scores
        lig_dir: str, default="ligand", subfolder to save extracted ligands
        molecule_dir: str,  default="molecule", subfolder to save extracted molecules (organic molecule, metal ion)
    """
    download_dir: str
    ccd_dir: str = "ccd"   
    ccd_path: str = "ccd/ccd.csv"  
    bcif_dir: str = "bcif"  
    cif_dir: str = "cif"   
    vs_dir: str = "entry"   
    lig_dir: str = "ligand"        
    molecule_dir: str = "molecule" 
    
    def __post_init__(self) -> None:
        for field in (fields(self)):
            if field.name != "download_dir":
                abs_path = os.path.join(self.download_dir, field.default)
                self.__setattr__(field.name, abs_path)
        
        
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
    data_preprocessor = DataPreprocessor(pdbid_path=args.pdbid_path, **download_config)
    pdb_ccd_instance_map = data_preprocessor.run()
    dataset_generator = DatasetGenerator(mode=args.mode,
                                         pdb_ccd_instance_map=pdb_ccd_instance_map,
                                         output_dir=output_dir,
                                         mmseqs_exec=args.mmseqs_exec,
                                         **download_config)
    dataset_generator.run()
