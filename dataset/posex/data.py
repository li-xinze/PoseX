import os
import json
import pymol
import shutil
import inspect
import itertools
import subprocess
import pandas as pd
import networkx as nx
from Bio import SeqIO
from rdkit import Chem
from rdkit import RDLogger
from functools import wraps
from functools import partial
from pymol import cmd, stored
from Bio.PDB.Model import Model
from multiprocessing import Pool
from collections import defaultdict
from posex import ccd
from rdkit.Geometry import Point3D
from posex.align import CrossAlignment
from Bio.SeqRecord import SeqRecord, Seq
from networkx.algorithms import bipartite
from typing import Callable, Literal, Any
from jinja2 import Environment, FileSystemLoader
RDLogger.DisableLog('rdApp.*')
from posex.utils import my_tqdm
from posex.utils import get_covalent_ligands, get_ligand_atoms, get_molecule_atoms, get_protein_atoms, \
    has_unknown_atoms, filter_ligand_with_distance, check_ETKDG, save_single_conformation, cif2seq, generate_conformation
from posex.utils import supress_stdout, run_in_tmp_dir
from posex.mmcif import cif_to_seq


PDB_BLOCKLIST = ["7X48", "7UYC", "7WJD", "7DB4", "6ZYU", "7W2W", "7ZSQ", "8AVA"]
CCD_BLOCKLIST = ["I8P", "5A3", "U71", "UEV"]
VALID_ELEMENTS = set(("H", "C", "O", "N", "P", "S", "F", "Cl"))
NUM_CPUS = 100



class DatasetGenerator():
    """
    A class for generating dataset

    Args:
        mode(Literal["cross_dock", "self_dock"]): dataset mode
        pdbid_path(str): path to the downloaded txt file contaning Entry IDs
        output_dir(str): folder to save the dataset
        mmseqs_exec(str): path to the MMseqs2
    """
    def __init__(self, 
                 mode: Literal["cross_dock", "self_dock"], 
                 pdb_ccd_instance_map: dict[str, dict[str, list[str]]],
                 output_dir: str, 
                 mmseqs_exec: str, 
                 **kwargs: Any) -> None:
        modes = ["cross_dock", "self_dock"]
        if mode not in modes:
            raise ValueError(f"mode={mode} not in {modes}")
        self.mode = mode
        vs_dir = kwargs.get("vs_dir")
        ccd_path = kwargs.get("ccd_path")
        components_path = kwargs.get("components_path")
        self.cif_dir = kwargs.get("cif_dir")
        self.lig_dir = kwargs.get("lig_dir")
        self.molecule_dir = kwargs.get("molecule_dir")
        self.output_dir = output_dir
        self.mmseqs_exec = mmseqs_exec
        self.pdb_ccd_instance_map = pdb_ccd_instance_map
        self.ccd_inchi_map, self.ccd_molwt_map = self._create_ccd_map(ccd_path)
        self.pdb_vs_map = self._create_pdb_vs_map(vs_dir)
        self.pdb_ccd_dict = self._create_pdb_ccd_dict()
        self.ccd_mol_dict = self._create_ccd_mol_dict()
        self.filter_record = self._init_filter_record()
        ccd.COMPONENTS_FILE = components_path
        

    def _init_filter_record(self):
        """ 
        """
        pdb_count = len(self.pdb_ccd_dict)
        ccd_count = len(set().union(*self.pdb_ccd_dict.values()))
        filter_record = {"input": {"PDB": pdb_count, "CCD": ccd_count}}
        print(f"PDB: {pdb_count}, CCD: {ccd_count}")
        return filter_record

    def _create_ccd_map(self, ccd_path: str) -> tuple[dict[str, str], dict[str, float]]:
        """Create a mapping from ccd to InChI and a mapping from ccd to MOLWT

        Args:
            ccd_path (str): path to save ccd table

        Returns:
            tuple[dict[str, str], dict[str, float]]: 
                a mapping from ccd to InChI and a mapping from ccd to MOLWT
        """
        df = pd.read_csv(ccd_path)
        df.CCD = df.CCD.apply(lambda x: x.strip("''"))
        df.set_index("CCD", drop=True, inplace=True)
        ccd_inchi_map = df.to_dict()["InChI"]
        ccd_molwt_map = df.to_dict()["MOLWT"]
        return ccd_inchi_map, ccd_molwt_map

    def _get_func_name(self, idx=1) -> str:
        """return the name of a funciton in the current call stack
        """
        return inspect.stack()[idx].function
    
    def _create_pdb_vs_map(self, vs_dir: str) -> dict[str, dict]:
        """Create a mapping from pdbid to validation info

        Args:
            pdbid_path (str): path to the downloaded txt file contaning Entry IDs
            vs_dir (str): folder to save json files containing validation scores

        Returns:
            dict[str, dict]: dict of pdbid to validation info
        """
        pdb_vs_map = defaultdict(set)
        for pdbid in self.pdb_ccd_instance_map:
            json_path = os.path.join(vs_dir, f"{pdbid}.json")
            with open(json_path, "r") as f:
                info = json.load(f)
                pdb_vs_map[pdbid] = info
        return pdb_vs_map

    def _create_pdb_ccd_dict(self) -> dict[str, set[str]]:
        """Create a dict of pdbid to ccd set
        """
        pdb_ccd_dict = defaultdict(set)
        for pdbid, ccd_info in self.pdb_ccd_instance_map.items():
            for ccd in ccd_info:
                pdb_ccd_dict[pdbid].add(ccd)
        return pdb_ccd_dict
            
    def _create_ccd_mol_dict(self) -> dict[str, Chem.Mol]:
        """Create a dict of ccd to rdkit mol
        """
        ccd_mol_dict = {}
        invalid_ccds = []
        ccd_set = set().union(*self.pdb_ccd_dict.values())
        for ccd in ccd_set:
            if ccd in self.ccd_inchi_map:
                ccd_mol_dict[ccd] = Chem.MolFromInchi(self.ccd_inchi_map[ccd], sanitize=False)
            else:
                invalid_ccds.append(ccd)
        # query
        print(f"invalid_ccds: {invalid_ccds}")
        return ccd_mol_dict

    def sync_dict(self) -> None:
        """ Synchronizing self.pdb_ccd_dict, self.ccd_mol_dict
        """
        ccd_keys = set(self.ccd_mol_dict.keys())
        self.pdb_ccd_dict = {k: v & ccd_keys for k, v in self.pdb_ccd_dict.items() if v & ccd_keys}
        ccd_values = set().union(*self.pdb_ccd_dict.values())
        self.ccd_mol_dict = {k: v for k, v in self.ccd_mol_dict.items() if k in ccd_values}
        
    def sync_filtered_result(func) -> Callable:
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            result = func(self, *args, **kwargs)
            self.sync_dict()
            func_name = func.__name__
            pdb_count = len(self.pdb_ccd_dict)
            ccd_count = len(self.ccd_mol_dict)
            print(f"{func_name:-<{50}} PDB: {pdb_count:-<{15}} CCD: {ccd_count:-<{15}}")
            self.filter_record[func_name] = {"PDB": pdb_count, "CCD": ccd_count}
            return result
        return wrapper
    
    @run_in_tmp_dir
    @sync_filtered_result
    @supress_stdout
    def filter_with_seq_length(self, max_len: int = 2000, num_cpus: int = NUM_CPUS) -> None:
        """filter with protein sequence length

        Args:
            max_len (int): the maximum acceptable length of protein sequence
        """
        pdbid_list = list(self.pdb_ccd_dict.keys())
        cif_to_seq_ = partial(cif_to_seq, cif_dir=self.cif_dir)
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(cif_to_seq_, pdbid_list)
            for pdbid, seq_dict in my_tqdm(pool_iter, total=len(pdbid_list), desc=self._get_func_name()):
                seq_len = 0
                chain_count = 0
                for item in seq_dict["sequences"]:
                    chain_count += len(item["protein"]["id"])
                    seq_len += len(item["protein"]["sequence"]) * len(item["protein"]["id"])
                    if len(item["protein"]["sequence"]) < 50:
                        del self.pdb_ccd_dict[pdbid]
                        break
                if seq_len > max_len or chain_count > 20:
                    if pdbid in self.pdb_ccd_dict:
                        del self.pdb_ccd_dict[pdbid]

    @sync_filtered_result
    def filter_with_unknown_ccd(self) -> None:
        """Remove unknown ligands (e.g. UNX, UNL)
        """
        pass

    @sync_filtered_result
    def filter_with_molwt(self, min_molwt: int = 100, max_molwt: int = 900) -> None:
        """Ligands weighing from 100 Da to 900 Da
        """
        self.ccd_mol_dict = {k: v for k, v in self.ccd_mol_dict.items() if  min_molwt <= self.ccd_molwt_map[k] < max_molwt}
    
    @sync_filtered_result
    def filter_with_num_heavy_atom(self, min_num_heavy_atom: int = 3) -> None:
        """Ligands with at least 3 heavy atoms
        """
        self.ccd_mol_dict = {k: v for k, v in self.ccd_mol_dict.items() if  v.GetNumHeavyAtoms() >= min_num_heavy_atom}

    @sync_filtered_result
    def filter_with_mol_element(self, valid_elements: set = VALID_ELEMENTS) -> None:
        """Ligands containing only H, C, O, N, P, S, F, Cl atoms
        """
        self.ccd_mol_dict = {k: v for k, v in self.ccd_mol_dict.items() 
                             if not set([atom.GetSymbol() for atom in v.GetAtoms()]) - valid_elements}

    @sync_filtered_result
    def filter_with_unknown_atoms(self, num_cpus: int = NUM_CPUS) -> None:
        """Structures with no unknown atoms (e.g. element X)
        """
        pdbid_list = list(self.pdb_ccd_dict.keys())
        has_unknown_atoms_ = partial(has_unknown_atoms, cif_dir=self.cif_dir)
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(has_unknown_atoms_, pdbid_list)
            for pdbid, flag in my_tqdm(pool_iter, total=len(pdbid_list), desc=self._get_func_name()):
                if flag:
                    del self.pdb_ccd_dict[pdbid]

    def filter_with_validation_report(self, indicator: str, is_valid_func: Callable) -> None:
        """ filter ccd instances with validation scores (e.g. RSR, RSCC, etc.)

        Args: 
            indicator: metric name (RSR, RSCC, completeness, stereo_outliers, intermolecular_clashes)
            is_valid_func: A function that takes metric value as input and returns boolean,
                           an instnace is valid if the return value is True
        """
        filtered_pdb_ccd_dict = {}
        for pdb, ccds in self.pdb_ccd_dict.items():
            for nonpolymer_entity in self.pdb_vs_map[pdb]["data"]["entry"]["nonpolymer_entities"]:
                ccd = nonpolymer_entity["nonpolymer_comp"]["chem_comp"]["id"]
                if ccd not in ccds: continue
                is_ccd_valid = True
                for instance in nonpolymer_entity["nonpolymer_entity_instances"]:
                    validation_score = instance["rcsb_nonpolymer_instance_validation_score"]
                    if validation_score is None: 
                        is_ccd_valid = False
                        break
                    score = validation_score[0][indicator]
                    if score is None or not is_valid_func(score):
                        is_ccd_valid = False
                        break
                if not is_ccd_valid:
                    ccds.remove(ccd)
            if ccds:
                filtered_pdb_ccd_dict[pdb] = ccds
        self.pdb_ccd_dict = filtered_pdb_ccd_dict

    @sync_filtered_result
    def filter_with_RSR(self, max_RSR: float = 0.3) -> None:
        """Ligand real space R-factor is at most 0.3
        """
        self.filter_with_validation_report(indicator="RSR", is_valid_func=lambda x: x < max_RSR)
       
    @sync_filtered_result
    def filter_with_RSCC(self, min_RSCC: float = 0.85) -> None:
        """Ligand real space correlation coefficient is at least 0.85
        """
        self.filter_with_validation_report(indicator="RSCC", is_valid_func=lambda x: x >= min_RSCC)

    @sync_filtered_result
    def filter_with_model_completeness(self) -> None:
        """Ligand model completeness is 100%
        """
        self.filter_with_validation_report(indicator="completeness", is_valid_func=lambda x: x == 1)

    @sync_filtered_result
    def filter_with_ETKDG(self, num_cpus: int = NUM_CPUS) -> None:
        """Ligand starting conformation could be generated with ETKDGv3
        """
        inputs = [(ccd, mol) for ccd, mol in self.ccd_mol_dict.items()]
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(check_ETKDG, inputs)
            for ccd, is_valid in my_tqdm(pool_iter, total=len(inputs), desc=self._get_func_name()):
                if not is_valid:
                    del self.ccd_mol_dict[ccd]

    @sync_filtered_result
    def filter_with_rdkit(self) -> None:
        """All ligand SDF files can be loaded with RDKit and pass its sanitization
        """
        pdbid_list = list(self.pdb_ccd_dict.keys())
        for pdbid in pdbid_list:
            ccds = self.pdb_ccd_dict[pdbid]
            filterd_ccds = ccds.copy()
            for ccd in ccds:
                ligand_path = f"{self.lig_dir}/{pdbid}_{ccd}.sdf"
                assert os.path.exists(ligand_path), ligand_path
                try:
                    sup = Chem.SDMolSupplier(ligand_path)
                    for mol in sup:
                        if mol is None:
                            filterd_ccds.remove(ccd)
                            break
                except Exception as e:
                    # print(e)
                    filterd_ccds.remove(ccd)
            if filterd_ccds:
                self.pdb_ccd_dict[pdbid] = filterd_ccds
            else:
                del self.pdb_ccd_dict[pdbid]

    @sync_filtered_result
    def filter_with_stereo_outliers(self) -> None:
        """PDB ligand report does not list stereochemical errors
        """
        self.filter_with_validation_report(indicator="stereo_outliers", is_valid_func=lambda x: x == 0)

    @sync_filtered_result
    def filter_with_intermolecular_clashes(self) -> None:
        """PDB ligand report does not list any atomic clashes
        """
        self.filter_with_validation_report(indicator="intermolecular_clashes", is_valid_func=lambda x: x == 0)

    @sync_filtered_result
    def filter_with_pdb_blocklist(self, blocklist: list = PDB_BLOCKLIST) -> None:
        """Blocklist for PDB entries
        """
        self.pdb_ccd_dict = {pdb: ccd for pdb, ccd in self.pdb_ccd_dict.items() if pdb not in blocklist}
        
    @sync_filtered_result
    def filter_with_ccd_blocklist(self, blocklist: list = CCD_BLOCKLIST) -> None:
        """Blocklist for CCD entries
        """
        self.ccd_mol_dict = {ccd: mol for ccd, mol in self.ccd_mol_dict.items() if ccd not in blocklist}
        
    @sync_filtered_result
    def filter_with_covalent_bond(self, num_cpus: int = NUM_CPUS) -> None:
        """Ligands that are not covalently bound to protein
        """
        get_covalent_ligands_ = partial(get_covalent_ligands, cif_dir=self.cif_dir)
        inputs = list(self.pdb_ccd_dict.items())
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(get_covalent_ligands_, inputs)
            for pdbid, covalent_ligands in my_tqdm(pool_iter, total=len(inputs), desc=self._get_func_name()):
                ccds = self.pdb_ccd_dict[pdbid]
                for ccd, asym_ids in self.pdb_ccd_instance_map[pdbid].items():
                    if ccd not in ccds: 
                        continue
                    for asym_id in asym_ids:
                        if (ccd, asym_id) in covalent_ligands:
                            ccds.remove(ccd)
                            break
                if ccds:
                    self.pdb_ccd_dict[pdbid] = ccds
                else:
                    del self.pdb_ccd_dict[pdbid]
                          
    @sync_filtered_result
    def filter_with_unique_pdb_ccd(self) -> None:
        """Get a set with unique pdbs and ccds by Hopcroft–Karp matching algorithm
        """
        G = nx.Graph()
        for pdb, ccd_set in self.pdb_ccd_dict.items():
            for ccd in ccd_set:
                G.add_edge(pdb, ccd)
        top_nodes = set(self.ccd_mol_dict.keys())
        matching = bipartite.maximum_matching(G, top_nodes)
        pdb_ccd_dict = {}
        for ccd in top_nodes:
            if ccd in matching:
                pdb_ccd_dict[matching[ccd]] = {ccd}
        self.pdb_ccd_dict = pdb_ccd_dict

    def filter_with_distance(self, 
                             min_dist: float, 
                             get_ligand_atom_func: Callable[[str, str], list] | None = None, 
                             get_protein_atom_func: Callable[[Model], list] | None = None, 
                             num_cpus: int = NUM_CPUS
                             ) -> None:
        """Filter out ligands that are not at the desired distance from other objects (ligand/protein)

        Args:
            min_dist (float): the minimum acceptable distance between a ligand and other objects
            get_ligand_atom_func (Callable[[str, str], list] | None): func to get ligand atoms
            get_protein_atom_func (Callable[[Model], list] | None)): func to get protein atoms
            num_cpus (int): num of processes
        """
        filter_ligand_with_distance_ = partial(filter_ligand_with_distance, 
                                               cif_dir=self.cif_dir, 
                                               min_dist=min_dist,
                                               get_ligand_atom_func=get_ligand_atom_func,
                                               get_protein_atom_func=get_protein_atom_func)
        inputs = list(self.pdb_ccd_dict.items())
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(filter_ligand_with_distance_, inputs)
            for pdbid, valid_ccds in my_tqdm(pool_iter, total=len(inputs), desc=self._get_func_name(2)):
                self.pdb_ccd_dict[pdbid] = valid_ccds

    @sync_filtered_result
    def filter_with_ligand_protein_distance(self, num_cpus: int = NUM_CPUS) -> None:
        """Intermolecular distance between the ligand(s) of interest and the protein is at least 0.2 Å
        """
        self.filter_with_distance(get_protein_atom_func=get_protein_atoms, min_dist=0.2, num_cpus=num_cpus)

    @sync_filtered_result
    def filter_with_ligand_ligand_distance(self, num_cpus: int = NUM_CPUS) -> None:
        """Intermolecular distance between ligand and other ligands is at least 5 Å
        """
        get_ligand_atoms_ = partial(get_ligand_atoms, cif_dir=self.cif_dir)
        self.filter_with_distance(get_ligand_atom_func=get_ligand_atoms_, min_dist=5, num_cpus=num_cpus)

    @sync_filtered_result
    def filter_with_ligand_organic_molecule_distance(self, num_cpus: int = NUM_CPUS) -> None:
        """Intermolecular distance between ligand(s) of interest and other small organic molecules is at least 0.2 Å
        """
        get_organic_molecule_atoms = partial(get_molecule_atoms, 
                                             molecule_dir=self.molecule_dir, 
                                             molecule_type="organic_molecule")
        self.filter_with_distance(get_ligand_atom_func=get_organic_molecule_atoms, min_dist=0.2, num_cpus=num_cpus)

    @sync_filtered_result
    def filter_with_ligand_metal_ion_distance(self, num_cpus: int = NUM_CPUS) -> None:
        """Intermolecular distance between the ligand(s) of interest and ion metals in complex is at least 0.2 Å
        """
        get_metal_ion_atoms = partial(get_molecule_atoms,
                                      molecule_dir=self.molecule_dir,
                                      molecule_type="metal_ion")
        self.filter_with_distance(get_ligand_atom_func=get_metal_ion_atoms, min_dist=0.2, num_cpus=num_cpus)

    @sync_filtered_result
    def select_single_conformation(self, num_cpus: int = NUM_CPUS) -> None:
        """Select single protein-ligand conformation
        """
        # for pdb entries
        single_conf_cif_dir = os.path.join(self.cif_dir, "single_conf")
        os.makedirs(single_conf_cif_dir, exist_ok=True)
        pdbid_list = list(self.pdb_ccd_dict.keys())
        save_single_conformation_ = partial(save_single_conformation, 
                                            cif_dir=self.cif_dir, 
                                            single_conf_cif_dir=single_conf_cif_dir)
        with Pool(processes=num_cpus) as pool:
            pool_iter = pool.imap_unordered(save_single_conformation_, pdbid_list)
            for _ in my_tqdm(pool_iter, total=len(pdbid_list), desc=self._get_func_name()):
                pass
        self.cif_dir = single_conf_cif_dir

    @run_in_tmp_dir
    @sync_filtered_result
    @supress_stdout
    def filter_with_crystal_contact(self, min_dist: int = 5) -> None:
        for pdbid, ccds in my_tqdm(self.pdb_ccd_dict.items(), desc=self._get_func_name()):
            cif_path = os.path.join(self.cif_dir, f"{pdbid}.cif")
            cmd.load(cif_path, pdbid)
            cmd.remove("solvent")
            valid_ccds = set()
            for ccd in ccds:
                cmd.select(ccd, f"not polymer and resn {ccd}")
                cmd.symexp("sym", pdbid, ccd, min_dist)
                objects = cmd.get_object_list()
                if len(objects) == 1:
                    valid_ccds.add(ccd)
                cmd.delete(ccd)
                cmd.delete("sym*")
            cmd.delete("all")
            self.pdb_ccd_dict[pdbid] = valid_ccds

    @run_in_tmp_dir
    @sync_filtered_result
    @supress_stdout
    def filter_with_clustering(self, num_cpus: int = NUM_CPUS) -> dict[str, set[str]]:
        """Cluster protein sequences

        Returns:
            dict[str, set[str]]: dict of the representative pdb to other pdbs in the cluster
        """
        res_name = "mmseqs_output"
        res_path = f"{res_name}_cluster.tsv"
        input_fasta_path = "input.fasta"
        pdbid_list = list(self.pdb_ccd_dict.keys())
        cif2seq_ = partial(cif2seq, cif_dir=self.cif_dir)
        
        with open(input_fasta_path, "w") as f:
            with Pool(processes=num_cpus) as pool:
                pool_iter = pool.imap_unordered(cif2seq_, pdbid_list)
                for pdbid, sequence in my_tqdm(pool_iter, total=len(pdbid_list), desc=self._get_func_name()):
                    record = SeqRecord(Seq(sequence), id=pdbid)
                    SeqIO.write(record, f, "fasta")

        mmseq2_script = f"{self.mmseqs_exec} easy-cluster {input_fasta_path} {res_name} tmp " + \
                        r"""--min-seq-id {min_seq_id} -c {c} --cluster-mode {cluster_mode}"""
        if self.mode == "cross_dock":
            running_script = mmseq2_script.format(
                min_seq_id=0.9,
                c=0.8,
                cluster_mode=1,
            )
        elif self.mode == "self_dock":
            running_script = mmseq2_script.format(
                min_seq_id=.0,
                c=1.0,
                cluster_mode=1,
            )
        else:
            raise ValueError()
        
        subprocess.run(running_script, shell=True, check=True, stdout=subprocess.DEVNULL)
        df = pd.read_csv(res_path, sep="\t", header=None)
        cluster_map = defaultdict(set)
        for _, row in df.iterrows():
            cluster_map[row[0]].add(row[1])
        if self.mode == "cross_dock":
            cluster_map = {k: v for k, v in cluster_map.items() if len(v) > 1}
            valid_pdbs = set().union(*cluster_map.values())
        else:
            pass
        # reselect representative_pdb
        new_cluster_map_except_k = {}
        for pdbs in cluster_map.values():
            pdbs = sorted(pdbs)
            representative_pdb = pdbs[0]
            new_cluster_map_except_k[representative_pdb] = pdbs[1:]
        if self.mode == "self_dock":
            self.pdb_ccd_dict = {pdb: ccd for pdb, ccd in self.pdb_ccd_dict.items() if pdb in new_cluster_map_except_k.keys()}
        elif self.mode == "cross_dock":
            self.pdb_ccd_dict = {pdb: ccd for pdb, ccd in self.pdb_ccd_dict.items() if pdb in valid_pdbs}
        else:
            raise ValueError()
        return new_cluster_map_except_k
        
    @run_in_tmp_dir
    @sync_filtered_result
    @supress_stdout
    def filter_with_cross_alignment(self, cluster_map: dict[str, set[str]]) -> set[tuple[tuple[str, str, str]]]:
        """Return a set of cross groups, each cross group contains tuple of (pdb, ccd, asym_id)

        Args: 
            cluster_map (dict[str, set[str]]): dict of the representative pdb to other pdbs in the cluster

        Returns:
            set[tuple[tuple[str, str, str]]]: a set of cross groups
        """
        cross_alignment = CrossAlignment(self.cif_dir, cluster_map, self.pdb_ccd_dict, self.pdb_ccd_instance_map)
        cross_groups = cross_alignment.run()
        pdb_ccd_dict = defaultdict(set)
        for cross_group in cross_groups:
            for pdb, ccd, ccd_asym_id in cross_group:
                pdb_ccd_dict[pdb].add(ccd)
        self.pdb_ccd_dict = pdb_ccd_dict
        return cross_groups

    def render_latex(self) -> None:
        output_path = f"{self.mode}_table.txt"
        env = Environment(loader=FileSystemLoader('template'))
        template = env.get_template(f"{self.mode}.txt")
        latex_output = template.render(**self.filter_record)
        with open(output_path, "w") as f:
            f.write(latex_output)

    @run_in_tmp_dir
    @supress_stdout
    def save_cross_dock_res(self, cross_groups: set[tuple[tuple[str, str, str]]]) -> None:
        """ Save cross-dock dataset

        Args: 
            cross_groups (set[tuple[tuple[str, str, str]]]):
                a set of cross groups, each cross group contains tuple of (pdb, ccd, asym_id)
        """
        os.makedirs(self.output_dir)
        for group_id, cross_group in  enumerate(cross_groups):
            pdb_num = len(cross_group)
            group_dir = os.path.join(self.output_dir, f"group_{group_id}_{pdb_num}")
            os.makedirs(group_dir)
            ref_pdb, _, _ = cross_group[0]
            ref_cif_path = os.path.join(self.cif_dir, f"{ref_pdb}.cif")
            for cand_pdb, cand_ccd, cand_ccd_asym_id in cross_group:
                ref_name = f"ref_prot_{cand_pdb}"
                cmd.load(ref_cif_path, ref_name)
                item_name = f"{cand_pdb}_{cand_ccd}"
                item_dir = os.path.join(group_dir, item_name)
                os.makedirs(item_dir)
                cif_path = os.path.join(self.cif_dir, f"{cand_pdb}.cif")
                cmd.load(cif_path, cand_pdb)
                cmd.align(f"{cand_pdb}////CA", f"{ref_name}////CA")
                cmd.remove("solvent")
                cmd.select(cand_ccd, f"not polymer and resn {cand_ccd} and {cand_pdb}")
                cmd.extract(cand_ccd, cand_ccd)
                cmd.delete(cand_ccd)
                ligand_sdf_path = f"{item_dir}/{item_name}_ligand.sdf"
                ligands_sdf_path = f"{item_dir}/{item_name}_ligands.sdf"
                raw_ligands_sdf_path = f"{self.lig_dir}/{item_name}.sdf"
                mols = list(Chem.SDMolSupplier(raw_ligands_sdf_path))
                cmd.load(raw_ligands_sdf_path, "ligands")
                # Transform ligands
                cmd.matrix_copy(cand_pdb, "ligands")
                stored.coords = []
                cmd.iterate_state(0, "ligands", "stored.coords.append((x, y, z))")
                i = 0
                for mol in mols:
                    conformer = mol.GetConformer()
                    for j in range(mol.GetNumAtoms()):
                        x, y, z = stored.coords[i]
                        i += 1
                        conformer.SetAtomPosition(j, Point3D(x, y, z))
                with Chem.SDWriter(ligands_sdf_path) as f:
                    for mol in mols:
                        f.write(mol)
                instance_name = f"{cand_pdb}_{cand_ccd}_{cand_ccd_asym_id}"
                selected_mol = None
                for mol in mols:
                    if mol.GetProp('_Name') == instance_name:
                        selected_mol = mol
                        break
                assert selected_mol is not None
                with Chem.SDWriter(ligand_sdf_path) as f:
                    f.write(selected_mol)

                pdb_save_path = f"{item_dir}/{item_name}_protein.pdb"
                cif_save_path = f"{item_dir}/{item_name}_protein.cif"
                cmd.save(pdb_save_path, cand_pdb)
                shutil.copy(cif_path, cif_save_path)
                cmd.delete("all")
        self._postprocess_cross_dock_res()

    def _postprocess_cross_dock_res(self) -> None:
        """Unfold cross-dock dataset
        """
        cross_dock_dir = self.output_dir
        unfold_cross_dock_dir = self.output_dir
        for group_name in my_tqdm(os.listdir(cross_dock_dir), desc=self._get_func_name()):
            group_dir = os.path.join(cross_dock_dir, group_name)
            pair_perms = list(itertools.permutations(os.listdir(group_dir), 2))
            for pair in pair_perms:
                ref_item, item = pair
                item_dir = os.path.join(group_dir, item)
                ref_item_dir = os.path.join(group_dir, ref_item)
                ref_ccd = ref_item.split("_")[1]
                pdb = item.split("_")[0]
                pair_name = f"{pdb}_{ref_ccd}"
                pair_dir = os.path.join(unfold_cross_dock_dir, pair_name)
                os.makedirs(pair_dir, exist_ok=True)
                ref_protein_pdb = os.path.join(ref_item_dir, f"{ref_item}_protein.pdb")
                ref_ligand_sdf = os.path.join(ref_item_dir, f"{ref_item}_ligand.sdf")
                ref_ligands_sdf = os.path.join(ref_item_dir, f"{ref_item}_ligands.sdf")
                protein_pdb = os.path.join(item_dir, f"{item}_protein.pdb")
                pair_ref_protein_pdb = os.path.join(pair_dir, f"{pair_name}_ref_protein.pdb")
                pair_protein_pdb = os.path.join(pair_dir, f"{pair_name}_protein.pdb")
                pair_ligand_sdf = os.path.join(pair_dir, f"{pair_name}_ligand.sdf")
                pair_ligands_sdf = os.path.join(pair_dir, f"{pair_name}_ligands.sdf")
                shutil.copy(ref_protein_pdb, pair_ref_protein_pdb)
                shutil.copy(protein_pdb, pair_protein_pdb)
                shutil.copy(ref_ligand_sdf, pair_ligand_sdf)
                shutil.copy(ref_ligands_sdf, pair_ligands_sdf)
                # write group id
                group_id_path = os.path.join(pair_dir, "group_id.txt")
                with open(group_id_path, "w") as f:
                    f.write(ref_item + "\n")
                    f.write(group_name)
                # write mol start conf 
                mol = Chem.SDMolSupplier(ref_ligand_sdf)[0]
                smiles = Chem.MolToSmiles(mol)
                ligand_start_conf_sdf_path = os.path.join(pair_dir, f"{pair_name}_ligand_start_conf.sdf")
                with Chem.SDWriter(ligand_start_conf_sdf_path) as f:
                    res, start_mol = generate_conformation(mol)
                    assert res != -1
                    f.write(start_mol)

                # save AF3 input json file
                _, task_dict = cif_to_seq(f"{item}_protein", item_dir)
                task_dict["sequences"].append({
                    "ligand": {
                        "id": "ZZ",
                        "smiles": smiles
                    }
                })
                task_dict["name"] = pair_name
                task_dict["modelSeeds"] = [42]
                task_dict["dialect"] = "alphafold3"
                task_dict["version"] = 1
                with open(os.path.join(pair_dir, f"{pair_name}.json"), "w") as f:
                    json.dump(task_dict, f, indent=2)
            shutil.rmtree(group_dir)
       
    @run_in_tmp_dir
    @supress_stdout
    def save_self_dock_res(self) -> None:
        """ Save self-dock dataset
        """
        cmd.reinitialize()
        os.makedirs(self.output_dir)
        for pdb, ccd_set in self.pdb_ccd_dict.items():
            for ccd in ccd_set:
                item_name = f"{pdb}_{ccd}"
                item_dir = os.path.join(self.output_dir, item_name)
                os.makedirs(item_dir)
                cif_path = os.path.join(self.cif_dir, f"{pdb}.cif")
                cmd.load(cif_path, pdb)
                cmd.remove("solvent")
                cmd.select(ccd, f"not polymer and resn {ccd} and {pdb}")
                cmd.extract(ccd, ccd)
                cmd.delete(ccd)
                ligand_sdf_path = f"{item_dir}/{item_name}_ligand.sdf"
                ligands_sdf_path = f"{item_dir}/{item_name}_ligands.sdf"
                ligand_start_conf_sdf_path = f"{item_dir}/{item_name}_ligand_start_conf.sdf"
                raw_ligands_sdf_path = f"{self.lig_dir}/{item_name}.sdf"
                shutil.copy(raw_ligands_sdf_path, ligands_sdf_path)
                pdb_save_path = f"{item_dir}/{item_name}_protein.pdb"
                cmd.save(pdb_save_path, pdb)
                cmd.delete("all")
                mol = Chem.SDMolSupplier(ligands_sdf_path)[0]
                smiles = Chem.MolToSmiles(mol)
                with Chem.SDWriter(ligand_sdf_path) as f:
                    f.write(mol)
                # write mol start conf 
                with Chem.SDWriter(ligand_start_conf_sdf_path) as f:
                    res, start_mol = generate_conformation(mol)
                    assert res != -1
                    f.write(start_mol)

                # save AF3 input json file
                _, task_dict = cif_to_seq(pdb, self.cif_dir)
                task_dict["sequences"].append({
                    "ligand": {
                        "id": "Z",
                        "smiles": smiles
                    }
                })
                task_dict["name"] = item_name
                task_dict["modelSeeds"] = [42]
                task_dict["dialect"] = "alphafold3"
                task_dict["version"] = 1
                with open(os.path.join(item_dir, f"{item_name}.json"), "w") as f:
                    json.dump(task_dict, f, indent=2)


                    
    def set_pdb_ccd_dict(self, pdb_ccd_dict) -> None:
        self.pdb_ccd_dict = pdb_ccd_dict

    def run(self):
        self.filter_with_unknown_ccd()
        self.filter_with_seq_length()
        self.filter_with_molwt()
        self.filter_with_num_heavy_atom()
        self.filter_with_mol_element()
        self.filter_with_covalent_bond()
        self.filter_with_unknown_atoms()
        self.filter_with_RSR()
        self.filter_with_RSCC()
        self.filter_with_model_completeness()
        self.filter_with_ETKDG()
        self.filter_with_rdkit()
        self.filter_with_stereo_outliers()
        self.filter_with_intermolecular_clashes()
        self.select_single_conformation()
        self.filter_with_ligand_protein_distance()
        if self.mode == "cross_dock":
            self.filter_with_ligand_ligand_distance()
            self.filter_with_crystal_contact()
            cluster_map = self.filter_with_clustering()
            cross_groups = self.filter_with_cross_alignment(cluster_map)
            self.save_cross_dock_res(cross_groups)
        elif self.mode == "self_dock":
            self.filter_with_ligand_organic_molecule_distance()
            self.filter_with_ligand_metal_ion_distance()
            self.filter_with_crystal_contact()
            self.filter_with_unique_pdb_ccd()
            self.filter_with_clustering()
            self.save_self_dock_res()
        else:
            raise ValueError(self.mode)
        # self.render_latex()