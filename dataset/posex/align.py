import os
import inspect
from pymol import cmd
from posex.utils import my_tqdm

class CrossAlignment():
    """
    A class for creating cross groups for cross-docking

    Args:
        cif_dir (str) : folder to save pdb entries (cif format)
        cluster_map (dict[str, set[str]]):   
            dict of the representative pdb to other pdbs in the cluster
        pdb_ccd_dict (dict[str, set[str]]):
            dict of PDBID to CCD set
        pdb_ccd_instance_map (dict[str, dict[str, list[str]]]): 
            dict of pdb to dict of ccd to list of asym_ids
        downsample_num (int): maximun num of items in a cross group

    """
    def __init__(self, 
                 cif_dir: str, 
                 cluster_map: dict[str, set[str]], 
                 pdb_ccd_dict: dict[str, set[str]],
                 pdb_ccd_instance_map: dict[str, dict[str, list[str]]], 
                 downsample_num: int = 8):
        self.cif_dir = cif_dir
        self.cluster_map = cluster_map
        self.pdb_ccd_dict = pdb_ccd_dict
        self.pdb_ccd_instance_map = pdb_ccd_instance_map   
        self.downsample_num = downsample_num

    def _get_func_name(self, idx=1) -> str:
        """Return the name of a funciton in the current call stack
        """
        return inspect.stack()[idx].function
    
    def check_lig_from_candidate_to_reference(self, 
                                              ref_pdb: str, 
                                              ref_ccd: str, 
                                              asym_ids: list, 
                                              valid_pdb_set: set
                                              ) -> dict[tuple[str, str], str]:
        """Check if each instance of the reference ligand is within 4.0 Å of the candidate ligand

        Args:
            ref_pdb (str): PDBID
            ref_ccd (str): CCD
            asym_ids (str): list of asym_ids
            valid_pdb_set (set): set of valid PDBIDs

        Returns:
            dict[tuple[str, str], str]: candidate_items(dict[tuple[pdb, ccd], ccd_asym_id])
        """
        candidate_items = {}
        for cand_pdb in sorted(valid_pdb_set):
            for cand_ccd in self.pdb_ccd_dict[cand_pdb]:
                if cand_ccd == ref_ccd: continue
                cand_ccd_asym_ids = self.pdb_ccd_instance_map[cand_pdb][cand_ccd]
                hit_asym_ids = set()
                for cand_ccd_asym_id in cand_ccd_asym_ids:
                    for asym_id in asym_ids:
                        ref_lig = f"{ref_pdb}_{ref_ccd}_{asym_id}"
                        res_name = f"{ref_lig}_{cand_ccd}_{cand_ccd_asym_id}"
                        cmd.select(ref_lig, f"resn {ref_ccd} and segi {asym_id} and {ref_pdb}")
                        cmd.select(res_name, f"resn {cand_ccd}  and segi {cand_ccd_asym_id} and ({cand_pdb} within 4.0 of {ref_lig})")
                        if cmd.count_atoms(res_name) > 0:
                            hit_asym_ids.add(asym_id)
                            if asym_id == asym_ids[0]:
                                selelcted_cand_ccd_asym_id = cand_ccd_asym_id
                        cmd.delete(res_name)
                        cmd.delete(ref_lig)
                if len(hit_asym_ids) == len(asym_ids):
                    candidate_items[(cand_pdb, cand_ccd)] = selelcted_cand_ccd_asym_id
        return candidate_items

    def check_lig_from_reference_to_candidate(self, 
                                              ref_pdb: str, 
                                              ref_ccd: str, 
                                              candidate_items: dict[tuple[str, str], str]
                                              ) -> dict[tuple[str, str], str]:
        """Check if each instance of the candidate ligand is within 4.0 Å of the reference ligand

        Args:
            ref_pdb (str): PDBID
            ref_ccd (str): CCD
            candidate_items (dict[tuple[str, str], str]): 
                dict[tuple[pdb, ccd], ccd_asym_id]

        Returns:
            dict[tuple[str, str], str]: filtered candidate_items
        """
        for cand_item in list(candidate_items.keys()):
            cand_pdb, cand_ccd = cand_item
            for asym_id in self.pdb_ccd_instance_map[cand_pdb][cand_ccd]:
                cand_lig = f"{cand_pdb}_{cand_ccd}_{asym_id}"
                res_name = f"{cand_lig}_{ref_ccd}"
                cmd.select(cand_lig, f"resn {cand_ccd} and segi {asym_id} and {cand_pdb}")
                cmd.select(res_name, f"resn {ref_ccd} and ({ref_pdb} within 4.0 of {cand_lig})")
                cmd.delete(cand_lig)
                if cmd.count_atoms(res_name) == 0:
                    del candidate_items[cand_item]
                    cmd.delete(res_name)
                    break
                else:
                    cmd.delete(res_name)
        return candidate_items
    
    def select_cross_candidates(self, 
                                ref_pdb: str, 
                                ref_ccd: str, 
                                asym_ids: list[str], 
                                valid_pdb_set: set
                                ) -> tuple[tuple[str, str, str]] | None:
        """Given a reference protien and ligand, return a cross group containing tuple of (pdb, ccd, asym_id)

        Args:
            ref_pdb (str): PDBID
            ref_ccd (str): CCD
            asym_ids (str): list of asym_ids
            valid_pdb_set (set): set of valid PDBIDs

        Returns:
            tuple[tuple[str, str, str]] | None: 
                tuple[tuple[pdb, ccd, ccd_asym_id]]
        """
        selected_ref_ccd_asym_id = asym_ids[0]
        candidate_items = self.check_lig_from_candidate_to_reference(ref_pdb, ref_ccd, asym_ids, valid_pdb_set)
        candidate_items = self.check_lig_from_reference_to_candidate(ref_pdb, ref_ccd, candidate_items)
        if candidate_items:
            # remove repeated ccd in candidate_items
            new_candidate_items = {}
            cand_ccds = set()
            for cand_item, cand_ccd_asym_id in candidate_items.items():
                cand_pdb, cand_ccd = cand_item
                if cand_ccd not in cand_ccds:
                    new_candidate_items[cand_item] = cand_ccd_asym_id
                    cand_ccds.add(cand_ccd)
            cross_group = [(ref_pdb, ref_ccd, selected_ref_ccd_asym_id)]
            for cand_item, cand_ccd_asym_id in new_candidate_items.items():
                cross_group.append((*cand_item, cand_ccd_asym_id))
            cross_group = tuple(sorted(cross_group))
        else:
            cross_group = None 
        return cross_group
    
    def filter_pdb_with_alignment(self, ref_pdb: str, pdb_set: set, max_rmsd: float = 2.0) ->  set:
        """Remove PDB entry if the RMSD of the protein alignment is greater than 2.0 Å

        Args:
            ref_pdb (str): reference PDBID
            pdb_set (set): PDBIDs of proteins to be aligned
            max_rmsd (str): the maximum acceptable alignment rmsd
            valid_pdb_set (set): set of valid PDBIDs

        Returns:
            set: valid pdb set
        """
        ref_pdb_path = os.path.join(self.cif_dir, f"{ref_pdb}.cif")
        cmd.load(ref_pdb_path, ref_pdb)
        valid_pdb_set = set()
        for pdb in pdb_set:
            pdb_path = os.path.join(self.cif_dir, f"{pdb}.cif")
            cmd.load(pdb_path, pdb)
            align_res = cmd.align(f"{pdb}////CA", f"{ref_pdb}////CA")
            rmsd = align_res[0]
            if rmsd <= max_rmsd:
                valid_pdb_set.add(pdb)
            else:
                cmd.delete(pdb)
        return valid_pdb_set
    
    def run(self) -> set[tuple[tuple[str, str, str]]]:
        """Return a set of cross groups, each cross group contains tuple of (pdb, ccd, asym_id)
        """
        cross_groups = set()
        for ref_pdb, pdb_set in my_tqdm(self.cluster_map.items(), desc=self._get_func_name(2)):
            valid_pdb_set = self.filter_pdb_with_alignment(ref_pdb, pdb_set)
            for ref_ccd in self.pdb_ccd_dict[ref_pdb]:
                asym_ids = self.pdb_ccd_instance_map[ref_pdb][ref_ccd]
                cross_group = self.select_cross_candidates(ref_pdb, ref_ccd, asym_ids, valid_pdb_set)
                if cross_group is not None:
                    # downsample
                    if len(cross_group) > self.downsample_num:
                        cross_group = cross_group[:self.downsample_num]
                    cross_groups.add(cross_group)
            cmd.delete("all")
        return cross_groups