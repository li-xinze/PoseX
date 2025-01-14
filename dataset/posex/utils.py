import os
import shutil
import inspect
import tempfile
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from Bio.PDB import is_aa
from functools import partial
from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from rdkit.Chem import AllChem
from multiprocessing import Pool
from typing import Callable, Literal
from Bio.PDB import MMCIFParser, MMCIFIO
import biotite.structure.io.pdbx as pdbx
from pdbecif.mmcif_io import CifFileReader
from Bio.PDB import Selection, NeighborSearch
from rcsb.utils.io.MarshalUtil import MarshalUtil

NUM_CPUS = 100
PERIODIC_TABLE = Chem.GetPeriodicTable()
ELEMENT_TABLE = set([PERIODIC_TABLE.GetElementSymbol(atomic_num).upper() for atomic_num in range(1, 119)])


class my_tqdm(tqdm):
    """A custom tqdm class
    """
    def __init__(self, *args, **kwargs):
        kwargs['bar_format'] = "{l_bar}{bar:5}{r_bar}"
        kwargs['leave'] = False
        super(my_tqdm, self).__init__(*args, **kwargs)
    

def create_pdb_ccd_instance_map(pdbid_list: list, 
                                cif_dir: str, 
                                num_cpus: int = NUM_CPUS
                                ) -> dict[str, dict[str, list[str]]]:
    """Create a dict of pdb to a dict of ccd to list of asym_ids

    Args:
        pdbid_list (list): list of PDBID
        cif_dir (str): folder containing pdb entries (cif format)

    Returns:
        dict[str, dict[str, list[str]]]:
            dict of pdb to dict of ccd to list of asym_ids
    """
    pdb_ccd_instance_map = {}
    get_ccd_instance_map_ = partial(get_ccd_instance_map, cif_dir=cif_dir)
    func_name = inspect.stack()[0].function
    with Pool(processes=num_cpus) as pool:
        pool_iter = pool.imap_unordered(get_ccd_instance_map_, pdbid_list)
        for pdbid, ccd_instance_map in my_tqdm(pool_iter, total=len(pdbid_list), desc=func_name):
            pdb_ccd_instance_map[pdbid] = ccd_instance_map
    return pdb_ccd_instance_map
    
def get_ccd_instance_map(pdbid, cif_dir):
    """Create a dict of ccd to list of asym_ids

    Returns:
        dict[str, list[str]]: dict of ccd to list of asym_ids
    """
    def _search_ligand_code(cif_data, entity_id):
        for r in cif_data._pdbx_entity_nonpoly.search('entity_id', entity_id).values():
            comp_id = r['comp_id']
        return comp_id
    
    def _search_asym_id(cif_data, entity_id):
        asym_ids = [r['id'] for r in cif_data._struct_asym.search('entity_id', entity_id).values()]
        return asym_ids

    cif_path = os.path.join(cif_dir, f"{pdbid}.cif")
    cfr = CifFileReader()
    ccd_instance_map = {}
    cif_obj = cfr.read(cif_path, output = 'cif_wrapper', ignore = ['_atom_site'])
    cif_data = list(cif_obj.values())[0]
    for entity_id, entity_type in zip(cif_data['_entity']['id'], cif_data['_entity']['type']): 
        if entity_type == 'non-polymer':
            ligand_code = _search_ligand_code(cif_data, entity_id)
            asym_ids = _search_asym_id(cif_data, entity_id)
            ccd_instance_map[ligand_code] = asym_ids
    return pdbid, ccd_instance_map

def run_in_tmp_dir(func: Callable) -> Callable:
    """A wrapper for running a function in tmp dir
    """
    def wrapper(*args, **kwargs):
        tmp_dir = next(tempfile._get_candidate_names())
        os.makedirs(tmp_dir, exist_ok=True)
        os.chdir(tmp_dir)
        try:
            result = func(*args, **kwargs)
        except:
            raise
        finally:
            os.chdir("..")
            shutil.rmtree(tmp_dir)
        return result
    return wrapper
    
def supress_stdout(func: Callable) -> Callable:
    """A wrapper for muting a funciton
    """
    def wrapper(*args, **kwargs):
        with open(os.devnull, 'w') as devnull:
            old_stdout = os.dup(1)
            os.dup2(devnull.fileno(), 1)
            try:
                result = func(*args, **kwargs)
            finally:
                os.dup2(old_stdout, 1)  
        return result
    wrapper.__name__ = func.__name__
    return wrapper

def bcif2cif(bcif_file: str, bcif_dir: str, cif_dir: str) -> tuple[str, bool]:
    """Convert bcif file to cif file
    """
    marsha_util = MarshalUtil()
    bcif_path = os.path.join(bcif_dir, bcif_file)
    cif_file = bcif_file.replace(".bcif", ".cif")
    cif_path = os.path.join(cif_dir, cif_file)
    if os.path.exists(cif_path):
        return bcif_file, True
    data = marsha_util.doImport(bcif_path, fmt="bcif")
    success = marsha_util.doExport(cif_path, data, fmt="mmcif")
    if not data:
        success = False
    return bcif_file, success

def generate_conformation(mol: Chem.Mol, max_attempts=0) -> tuple[int, Chem.Mol]:
    """Generate conformation with ETKDGv3

    Args: 
        max_attempts (int): the maximum number of attempts to try embedding 

    Returns:
        tuple[int, Chem.Mol]: 
            tuple of new conformation ID and mol
    """
    mol.RemoveAllConformers()
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = False
    params.randomSeed = 42
    params.maxAttempts = max_attempts
    try:
        res = AllChem.EmbedMolecule(mol, params)
    except:
        res = -1
    return res, mol

def check_ETKDG(inputs, max_attempts: int = 0) -> tuple[str, bool]:
    """Check if starting conformation could be generated with ETKDGv3
    """
    ccd, input_mol = inputs
    input_mol = Chem.RemoveHs(input_mol)
    res, mol = generate_conformation(input_mol, max_attempts)
    mol = Chem.RemoveHs(mol)
    is_valid = (res != -1) & (input_mol.GetNumAtoms() == mol.GetNumAtoms())
    return ccd, is_valid
    
def has_unknown_atoms(pdbid, cif_dir) -> tuple[str ,bool]:
    """Check if structure has no unknown atoms
    """
    parser = MMCIFParser(QUIET=True)
    cif_path = os.path.join(cif_dir, f"{pdbid}.cif")
    protein = parser.get_structure("protein", cif_path)[0]
    atoms = Selection.unfold_entities(protein, 'A') 
    for atom in atoms:
        if atom.element not in ELEMENT_TABLE:
            return pdbid, True
    return pdbid, False

def cif2seq(pdbid, cif_dir, chain_sep="") -> tuple[str, str]:
    """Extract protein sequence from a cif file

    Args:
        pdbid (str): PDBID
        cif_dir (str): folder containing pdb entries (cif format)
        chain_sep (str): delimiter between different chains

    Returns:
        tuple[str, str]: PDBID, sequence
    """
    seq_dict = {}
    cif_path = f"{cif_dir}/{pdbid}.cif"
    with open(cif_path, "rt") as f:
        cif = pdbx.CIFFile.read(f)
    category = cif.block["entity_poly"]
    category_dict = {k: column.as_array() for k, column in category.items()}
    df_entity_poly = pd.DataFrame(category_dict, dtype=str)
    seqs = df_entity_poly["pdbx_seq_one_letter_code_can"].values
    entity_ids = df_entity_poly["entity_id"].values
    pdbx_strand_id = df_entity_poly["pdbx_strand_id"].values
    for entity_id, chains, seq in zip(entity_ids, pdbx_strand_id, seqs):
        for chain in chains.split(","):
            seq_dict[f"{entity_id}_{chain}"] = seq.replace("\n", "")
    seq = chain_sep.join(seq_dict.values())
    return pdbid, seq

def save_single_conformation(pdbid: str, cif_dir: str, single_conf_cif_dir: str, overwrite: bool = True) -> None:
    """Save the first alternative locations

    Args:
        pdbid (str): PDBID
        cif_dir (str): folder containing pdb entries (cif format)
        single_conf_cif_dir (str): folder to save processed pdb entries 
        overwrite (bool): if overwrite files
    """
    cif_path = os.path.join(cif_dir, f"{pdbid}.cif")
    single_conf_cif_path = os.path.join(single_conf_cif_dir, f"{pdbid}.cif")
    if not overwrite and os.path.exists(single_conf_cif_path):
        return
    parser = MMCIFParser(QUIET=True, auth_chains=False)
    try:
        protein = parser.get_structure("protein", cif_path)
    except KeyError:
        shutil.copy(cif_path, single_conf_cif_path)
        return 
    atoms = Selection.unfold_entities(protein, 'A') 
    valid_serial_numbers = []
    for atom in atoms:
        if atom.is_disordered() == 0:
            valid_serial_numbers.append(atom.serial_number)
        else:
            selected_disordered_id = atom.disordered_get_id_list()[0]
            valid_serial_numbers.append(atom.disordered_get(selected_disordered_id).serial_number)
    valid_serial_numbers = set(valid_serial_numbers)
    for key in list(parser._mmcif_dict.keys()):
        if key.startswith('_atom_site.'):
            new_values = [v for idx, v in enumerate(parser._mmcif_dict[key]) if idx + 1 in valid_serial_numbers]
            parser._mmcif_dict[key] = new_values
    io = MMCIFIO()
    io.set_dict(parser._mmcif_dict)
    io.save(single_conf_cif_path, preserve_atom_numbering=True)

def get_protein_atoms(protein: Model, ligands: set[tuple[str, str]]) -> list[Atom]:
    """Get protein atoms in a pdb entry

    Args:
        protein (Model): Bio.PDB.Model.Model
        ligands (tuple[str, str]): set of ligands in form of (ccd, asym_id)
       
    Returns:
        list[Atom]: atoms
    """
    residues = Selection.unfold_entities(protein, 'R') 
    atoms = []
    for residue in residues:
        ccd = residue.get_resname()
        asym_id = residue.get_full_id()[2]
        ligand = (ccd, asym_id)
        if is_aa(residue) and ligand not in ligands:
            atoms += list(residue.get_atoms())
    return atoms

def get_ligand_atoms(pdbid: str, except_ccd: str, cif_dir: str) -> list[Atom]:
    """Get ligand atoms in a pdb entry
    
    Args:
        pdbid (str): PDBID
        except_ccd (set): ccd of excluded ligands
        cif_dir (str): folder containing pdb entries (cif format)

    Returns:
        list[Atom]: atoms
    """
    cif_path = os.path.join(cif_dir, f"{pdbid}.cif")
    parser = MMCIFParser(QUIET=True)
    protein = parser.get_structure("protein", cif_path)[0]
    residues = Selection.unfold_entities(protein, 'R') 
    except_ccds = (except_ccd, "HOH", "WAT")
    ligand_atoms = []
    for residue in residues:
        if residue.get_resname() not in except_ccds and residue.get_id()[0] != ' ':
            for atom in residue.get_atoms():
                ligand_atoms.append(atom)
    return ligand_atoms

def get_molecule_atoms(pdbid: str, 
                       except_ccd: set, 
                       molecule_dir: str, 
                       molecule_type: Literal["organic_molecule", "metal_ion"]
                       ) -> list[Atom]:
    """Get atoms of molecules (organic molecule, metal ion)

    Args:
        pdbid (str): PDBID
        except_ccd (set): ccd of excluded ligands
        molecule_dir (str): folder containing extracted molecules (organic molecule, metal ion)
        molecule_type (Literal["organic_molecule", "metal_ion"]))

    Returns:
        list[Atom]: atoms
    """
    molecule_path = os.path.join(molecule_dir, f"{pdbid}_{molecule_type}.cif")
    parser = MMCIFParser(QUIET=True)
    try:
        protein = parser.get_structure("protein", molecule_path)[0]
    except KeyError:
        return []
    atoms = []
    for chain in protein.get_chains():
        for res in chain.get_residues():
            if res.get_resname() != except_ccd:
                for atom in res.get_atoms(): 
                    atoms.append(atom)
    return atoms

def filter_ligand_with_distance(inputs: tuple[str, set], 
                                cif_dir: str,
                                min_dist: float, 
                                get_ligand_atom_func: Callable[[str, str], list] | None, 
                                get_protein_atom_func: Callable[[Model], list] | None
                                ) -> set:
    """Filter out ligands that are not at the desired distance from other objects (ligand/protein)

    Args:
        inputs (tuple[str, set]): pdbid, ccd_set
        cif_dir (str): folder to save pdb entries (cif format)
        min_dist (float): the minimum acceptable distance between a ligand and other objects
        get_ligand_atom_func (Callable[[str, str], list] | None): func to get ligand atoms
        get_protein_atom_func (Callable[[Model], list] | None)): func to get protein atoms

    Returns:
        set: invalid ccd set
    """
    pdbid, ccds = inputs
    valid_ccds = set()
    cif_path = os.path.join(cif_dir, f"{pdbid}.cif")
    parser = MMCIFParser(auth_chains=False)
    protein = parser.get_structure("protein", cif_path)[0]
    if get_protein_atom_func is not None:
        _, ccd_instance_map = get_ccd_instance_map(pdbid, cif_dir)
        ligands = set((ccd, asym_id) for ccd, asym_ids in ccd_instance_map.items() 
                      for asym_id in asym_ids)
        atoms = get_protein_atom_func(protein, ligands)
        ns = NeighborSearch(atoms)
    for ccd in ccds:
        is_geq_limit = False
        for chain in protein.get_chains():
            for residue in chain.get_residues():
                if ccd == residue.get_resname():
                    if get_ligand_atom_func is not None:
                        atoms = get_ligand_atom_func(pdbid, ccd)
                        if not atoms: 
                            continue
                        ns = NeighborSearch(atoms)
                    for atom in residue.get_atoms(): 
                        if ns.search(atom.coord, min_dist, 'A'):
                            is_geq_limit = True
                            break
        if not is_geq_limit:
            valid_ccds.add(ccd)               
    return pdbid, valid_ccds

def get_covalent_ligands(inputs: tuple[str, set], cif_dir: str) -> tuple[str, list]:
    """Extract covalent ligands from a given ccd set

    Args:
        inputs (tuple[str, set]): pdbid, ccd_set
        cif_dir (str): folder to save pdb entries (cif format)

    Returns:
        tuple[str, list]: pdbid, list of tuple(ligand_ccd, asym_id)
    """
    pdbid, ccds = inputs
    cif_path = os.path.join(cif_dir, f"{pdbid}.cif")
    def search_asym_id(cif_data, entity_id):
        asym_id = []
        if '_struct_asym' in cif_data:
            for r in cif_data._struct_asym.search('entity_id',entity_id).values():
                asym_id.append(r['id'])
            asym_id = ','.join(asym_id)
        return asym_id
    
    cfr = CifFileReader()
    cif_obj = cfr.read(cif_path, output = 'cif_wrapper', ignore = ['_atom_site']) 
    cif_data = list(cif_obj.values())[0]
    polymer_chain_set = set()
    ligands_with_covalent_bond = []
    entity = cif_data['_entity']
    for entity_id, entity_type in zip(entity['id'], entity['type']):
        asym_id = search_asym_id(cif_data, entity_id)
        if entity_type == 'polymer': 
            polymer_chain_set.add(asym_id)

    if '_struct_conn' in cif_data:
        for r in cif_data._struct_conn.search('conn_type_id', 'covale').values():
            if (r.get("ptnr2_label_comp_id") in ccds and r.get('ptnr1_label_asym_id') in polymer_chain_set):
                ligand_ccd = r.get("ptnr2_label_comp_id")
                asym_id = r.get('ptnr2_label_asym_id')
                ligands_with_covalent_bond.append((ligand_ccd, asym_id))
            elif (r.get("ptnr1_label_comp_id") in ccds and r.get('ptnr2_label_asym_id') in polymer_chain_set):
                ligand_ccd = r.get("ptnr1_label_comp_id")
                asym_id = r.get('ptnr1_label_asym_id')
                ligands_with_covalent_bond.append((ligand_ccd, asym_id))
            else:
                pass
    return pdbid, ligands_with_covalent_bond