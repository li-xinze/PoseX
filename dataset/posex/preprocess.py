import os
import json
import inspect
import requests
import subprocess
import pandas as pd
from pymol import cmd
from typing import Any
from rdkit import Chem
from functools import partial
from multiprocessing import Pool
from collections import defaultdict
from posex.utils import bcif2cif, my_tqdm
from pdbecif.mmcif_io import CifFileReader
from posex.utils import run_in_tmp_dir, create_pdb_ccd_instance_map


NUM_CPUS = 100
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class DataPreprocessor():
    """
    A class for downloading and preprocessing data

    Args:
        pdbid_path : path to the downloaded txt file contaning Entry IDs
    """
    def __init__(self, pdbid_path: str, **kwargs: Any):
        self.bcif_dir = kwargs.get("bcif_dir")
        self.cif_dir = kwargs.get("cif_dir")
        self.vs_dir = kwargs.get("vs_dir")
        self.ccd_dir = kwargs.get("ccd_dir")
        self.ccd_path = kwargs.get("ccd_path")
        self.lig_dir = kwargs.get("lig_dir")
        self.molecule_dir = kwargs.get("molecule_dir")
        with open(pdbid_path, "r") as f:
            self.pdbid_list = f.readline().strip().split(",")
        

    def _get_func_name(self, idx=1) -> str:
        """return the name of a funciton in the current call stack
        """
        return inspect.stack()[idx].function

    def download_bcif(self) -> None:
        os.makedirs(self.bcif_dir, exist_ok=True)
        for pdbid in my_tqdm(self.pdbid_list, desc=self._get_func_name()):
            if os.path.exists(f"{self.bcif_dir}/{pdbid}.bcif"):
                continue
            try:
                url = f"https://models.rcsb.org/{pdbid}.bcif"
                subprocess.run(["wget", url], cwd=self.bcif_dir, check=True)
            except Exception as e:
                print(pdbid, e)

    def convert_bif_to_cif(self, num_cpus=NUM_CPUS) -> None:
        os.makedirs(self.cif_dir, exist_ok=True)
        bcif2cif_ = partial(bcif2cif, bcif_dir=self.bcif_dir, cif_dir=self.cif_dir)
        with Pool(processes=num_cpus) as pool:
            for bcif_file, success in my_tqdm(pool.imap_unordered(bcif2cif_, os.listdir(self.bcif_dir)), desc=self._get_func_name()):
                if not success:
                    print(f" {bcif_file} failed")
              
    def download_validation_score(self) -> None:
        os.makedirs(self.vs_dir, exist_ok=True)
        url = 'https://data.rcsb.org/graphql'
        with open(f"{PROJECT_DIR}/template/vs_query.txt", "r") as f:
            lines = f.readlines()
            query_string = "".join([line.strip() for line in lines])
        for pdbid in my_tqdm(self.pdbid_list, desc=self._get_func_name()):
            json_path = f"{self.vs_dir}/{pdbid}.json"
            if os.path.exists(json_path):
                continue
            content = {
                "query": query_string,
                "variables": {"id": pdbid}
            }
            response = requests.post(url, json=content)
            with open(json_path, 'w') as f:
                json.dump(response.json(), f)

    def query_ccd(self, ccd: str) -> tuple[float, str]:
        """Query formula_weight and inchi of a ligand using rcsb API

        Args:
            ccd (str): chemical descriptions of a ligand

        Returns:
            tuple[float, str]: formula_weight and inchi of the ligand
        """
        json_path = os.path.join(self.ccd_dir, f"{ccd}.json")
        if not os.path.exists(json_path):
            url = 'https://data.rcsb.org/graphql'
            with open(f"{PROJECT_DIR}/template/ccd_query.txt", "r") as f:
                lines = f.readlines()
            query_string = "".join([line.strip() for line in lines])
            content = {
                "query": query_string,
                "variables": {"id": ccd}
            }
            response = requests.post(url, json=content)
            with open(json_path, 'w') as f:
                json.dump(response.json(), f)

        with open(json_path, "r") as f:
            info = json.load(f)
            formula_weight = info["data"]["chem_comp"]["chem_comp"]["formula_weight"]
            inchi = None
            pdbx_chem_comp_descriptor = info["data"]["chem_comp"]["pdbx_chem_comp_descriptor"]
            if pdbx_chem_comp_descriptor is not None:
                for item in pdbx_chem_comp_descriptor:
                    if item["type"] == "InChI":
                        inchi = item["descriptor"]
                        break
        return formula_weight, inchi

    def _download_ccd_dict(self) -> str:
        """Download components.cif
        """
        os.makedirs(self.ccd_dir, exist_ok=True)
        components_path = os.path.join(self.ccd_dir, "components.cif")
        if not os.path.exists(components_path):
            print("downloading components.cif...")
            url = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif"
            subprocess.run(["wget", url], cwd=self.ccd_dir, check=True)
            print("finish")
        component_path = os.path.join(self.ccd_dir, "components.cif")
        print("reading components.cif...")
        print("finish")
        return component_path

    def build_ccd_table(self) -> None:
        """Create a ccd table and save it to self.ccd_path
        """
        if not os.path.exists(self.ccd_path):
            component_path = self._download_ccd_dict()
            cfr = CifFileReader()
            cif_obj = cfr.read(component_path, ignore = ['_atom_site'])
            df_dict = defaultdict(list)
            for ccd, values in cif_obj.items():
                if ccd == "UNL": continue
                descriptor_types = values["_pdbx_chem_comp_descriptor"]["type"]
                inchi_idx = -1
                for idx, dtype in enumerate(descriptor_types):
                    if dtype == "InChI":
                        inchi_idx = idx
                        break
                if inchi_idx == -1:
                    print(f"no InChI for {ccd}")
                else:
                    inchi = values["_pdbx_chem_comp_descriptor"]["descriptor"][inchi_idx]
                    df_dict["CCD"].append(ccd)
                    df_dict["InChI"].append(inchi)
                    df_dict["MOLWT"].append(values["_chem_comp"]["formula_weight"])
        else:
            df_dict = pd.read_csv(self.ccd_path).to_dict(orient="list")
            df_dict["CCD"] = [ccd.strip("''") for ccd in df_dict["CCD"]]
        
        # get all ccds
        ccd_set = set()
        for pdbid in self.pdbid_list:
            json_path = os.path.join(self.vs_dir, f"{pdbid}.json")
            with open(json_path, "r") as f:
                info = json.load(f)
                try:
                    for nonpolymer_entity in info["data"]["entry"]["nonpolymer_entities"]:
                        ccd = nonpolymer_entity["nonpolymer_comp"]["chem_comp"]["id"]
                        ccd_set.add(ccd)
                except:
                    print(json_path)
                    raise

        uncover_ccds = ccd_set - set(df_dict["CCD"])
        for ccd in my_tqdm(uncover_ccds, desc=self._get_func_name()):
            formula_weight, inchi = self.query_ccd(ccd)
            if formula_weight is not None and inchi is not None:
                df_dict["CCD"].append(ccd)
                df_dict["InChI"].append(inchi)
                df_dict["MOLWT"].append(formula_weight)
        df = pd.DataFrame(df_dict)
        
        df.CCD = df.CCD.apply(lambda x: f"\'{x}\'")
        df.to_csv(self.ccd_path, na_rep=None, index=False)

    @run_in_tmp_dir
    def extract_ligand(self) -> dict[str, dict[str, list[str]]]:
        """Extract sdf files of all ligands from cif files

        Returns:
            dict[str, dict[str, list[str]]]:
                dict of pdb to dict of ccd to list of asym_ids
        """
        pdb_ccd_instance_map = create_pdb_ccd_instance_map(self.pdbid_list, self.cif_dir)
        os.makedirs(self.lig_dir, exist_ok=True)
        for pdb, ccd_dict in my_tqdm(pdb_ccd_instance_map.items(), desc=self._get_func_name()):
            cif_path = os.path.join(self.cif_dir, f'{pdb}.cif')
            skip = True
            for ccd in ccd_dict:
                ligand_path = f"{self.lig_dir}/{pdb}_{ccd}.sdf"
                if not os.path.exists(ligand_path):
                    skip = False
                    break
            if skip:
                continue
            cmd.load(cif_path, pdb)
            for ccd, asym_ids in ccd_dict.items():
                ligand_path = f"{self.lig_dir}/{pdb}_{ccd}.sdf"
                mols = []
                for asym_id in asym_ids:
                    ligand_id = f"ligand_{asym_id}"
                    tmp_ligand_path = f"{self.lig_dir}/{ligand_id}.sdf"
                    cmd.select(ligand_id, f"resn {ccd} and segi {asym_id} and {pdb}")
                    cmd.save(tmp_ligand_path, ligand_id)
                    mol = Chem.SDMolSupplier(tmp_ligand_path)[0]
                    if mol is not None:
                        mols.append(mol)
                    else:
                        # All ligand SDF files can be loaded with RDKit and pass its sanitization
                        mols = []
                        break
                    os.remove(tmp_ligand_path)
                    cmd.delete(ligand_id)
                w = Chem.SDWriter(ligand_path)
                for mol in mols:
                    w.write(mol)
                w.close()
            cmd.delete("all")
        return pdb_ccd_instance_map

    @run_in_tmp_dir
    def extract_molecule(self) -> None:
        """Extract sdf files of organic_molecule and metal_ion from cif files
        """
        os.makedirs(self.molecule_dir, exist_ok=True)
        for pdb in my_tqdm(self.pdbid_list, desc=self._get_func_name()):
            cif_path = os.path.join(self.cif_dir, f'{pdb}.cif')
            organic_molecule_path = f"{self.molecule_dir}/{pdb}_organic_molecule.cif"
            metal_ion_path = f"{self.molecule_dir}/{pdb}_metal_ion.cif"
            if os.path.exists(organic_molecule_path) and os.path.exists(metal_ion_path):
                continue
            cmd.load(cif_path, pdb)
            cmd.select("organic_molecule", "organic")
            cmd.select("metal_ion", "metal")
            cmd.save(organic_molecule_path, "organic_molecule")
            cmd.save(metal_ion_path, "metal_ion")
            cmd.delete("all")

    def run(self) -> dict[str, dict[str, list[str]]]:
        """Preprocess

        Returns:
            dict[str, dict[str, list[str]]]:
                dict of pdb to dict of ccd to list of asym_ids
        """
        self.download_bcif()
        self.convert_bif_to_cif()
        self.download_validation_score()
        self.build_ccd_table()
        self.extract_molecule()
        pdb_ccd_instance_map = self.extract_ligand()
        return pdb_ccd_instance_map
