import sys

sys.path.append("./")
import numpy as np
from rdkit import Chem
from biotite.structure import AtomArray, concatenate
from biotite.structure.io import pdbx
from pathlib import Path
import json
import click
import pandas as pd
from tqdm import tqdm
from typing import *
from dataset.utils.repair_pdb import fix_pdb
from dataset.utils.openmm_helper import ProLigRelax
from dataset.utils.common_helper import create_logger
import traceback

logger = create_logger(__name__)
METALS = ["Na", "K", "Ca", "Mg", "Fe", "Zn", "Cu", "Mn", "Co", "Ni"]


def _format_4letter(atom_name: str):
    output = None
    if len(atom_name) == 4:
        output = atom_name
    elif len(atom_name) == 3:
        output = " " + atom_name
    elif len(atom_name) == 2:
        output = " " + atom_name + " "
    elif len(atom_name) == 1:
        output = " " + atom_name + "  "
    else:
        raise ValueError()

    return output


def _get_template_mol(cur_dir: Path):
    json_fn = list(cur_dir.glob("inputs.json"))[0]
    with open(json_fn, "r") as f:
        data = json.load(f)

    for sequence in data[0]["sequences"]:
        if "ligand" in sequence:
            smiles = sequence["ligand"]["ligand"]
            print(f"The first {smiles=}")
            break

    params = Chem.SmilesParserParams()
    lig_mol = Chem.MolFromSmiles(smiles, params)
    Chem.Kekulize(lig_mol, clearAromaticFlags=True)
    return lig_mol


def rdmol_to_file(mol: Chem.Mol, out_fn: Path | str, **kwargs) -> None:
    out_fn.parent.mkdir(exist_ok=True, parents=True)
    if out_fn.suffix == ".sdf":
        with Chem.SDWriter(out_fn) as writer:
            writer.write(mol)
    elif out_fn.suffix == ".mae":
        with Chem.MaeWriter(out_fn) as writer:
            writer.write(mol)
    elif out_fn.suffix == ".pdb":
        metals_uc = [x.upper() for x in METALS]
        with open(out_fn, "w") as w:
            w.write("REMARK 2025\n")
            w.write("REMARK 2025 WRITTEN BY Rdkit MolToPDBBlock\n")
            if kwargs.get("seq_blocks", None):
                [w.write(block) for chain, block in kwargs["seq_blocks"].items()]
            pdb_block = Chem.MolToPDBBlock(mol)
            pdb_lines = []
            prev_chain_name = None
            for idx, line in enumerate(pdb_block.split("\n")):
                if idx == 0 and prev_chain_name is None:
                    prev_chain_name = line[20:23].strip()

                if "<" in line:
                    line = line.replace("<", " ")

                if line[17:20].strip() in metals_uc + ["HOH"]:
                    line = "HETATM" + line[6:]

                cur_chain_name = line[20:22].strip()
                if cur_chain_name != prev_chain_name and line[:6] != "CONECT":
                    pdb_lines.append("TER")
                    prev_chain_name = cur_chain_name
                pdb_lines.append(line)

            w.write("\n".join(pdb_lines) + "\n")
    else:
        raise NotImplementedError()


def _calc_dist(coords1: np.ndarray, coords2: np.ndarray):
    dist = np.sqrt(np.sum((coords1[:, np.newaxis, :] - coords2) ** 2, axis=-1))
    return dist


def _get_charge(atom_name: str, res_name: str, label_charge: None) -> int:
    if atom_name == "OXT":
        charge = -1
    elif res_name == "LYS" and atom_name == "NZ":
        charge = 1
    elif res_name == "ARG" and atom_name == "NH2":
        charge = 1
    elif res_name == "ASP" and atom_name == "OD2":
        charge = -1
    elif res_name == "GLU" and atom_name == "OE2":
        charge = -1
    elif (label_charge is not None) and (str(label_charge) not in ["0", "?"]):
        charge = label_charge
    else:
        charge = 0
    return charge


def _add_conformer_to_mol(new_mol: Chem.Mol, pos: np.array) -> Chem.Mol:
    conf = Chem.Conformer(new_mol.GetNumAtoms())
    conf.SetPositions(pos)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conf, assignId=True)
    return new_mol


def extract_protein_and_ligand_array(this_array: AtomArray, lig_name: str) -> Chem.Mol:
    # prot_array: AtomArray, het_array: AtomArray, lig_array: AtomArray
    prot_array = this_array[this_array.hetero == False]
    het_array = this_array[this_array.hetero == True]
    lig_array = het_array[het_array.res_name == lig_name]
    lig_chain_id = np.unique(lig_array.chain_id)[0]
    return prot_array, lig_array, lig_chain_id


def _prot_atomarray_to_mol(prot_array: AtomArray) -> Chem.Mol:
    new_mol = Chem.RWMol()
    for idx, (
        element_name,
        res_name,
        res_id,
        chain_id,
        atom_name,
        ins_code,
        cif_charge,
    ) in enumerate(
        zip(
            prot_array.element,
            prot_array.res_name,
            prot_array.res_id,
            prot_array.chain_id,
            prot_array.atom_name,
            prot_array.ins_code,
            prot_array.charge,
        )
    ):
        if res_name == "HOH":
            pass
        rdatom = Chem.Atom(element_name)
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(f"{res_name:<4}")
        mi.SetResidueNumber(int(res_id))
        mi.SetChainId(str(chain_id[0]))
        mi.SetName(_format_4letter(atom_name))
        mi.SetInsertionCode("<>" if ins_code == "" else ins_code)
        # mi.SetInsertionCode('1')
        rdatom.SetMonomerInfo(mi)

        chg = _get_charge(atom_name, res_name, cif_charge)
        rdatom.SetFormalCharge(int(chg))
        atom_idx = new_mol.AddAtom(rdatom)

    for at1, at2, bt in prot_array.bonds.as_array():
        if bt == 1:
            new_mol.AddBond(int(at1), int(at2), Chem.BondType.SINGLE)
        elif bt == 2:
            new_mol.AddBond(int(at1), int(at2), Chem.BondType.DOUBLE)
        elif bt == 3:
            new_mol.AddBond(int(at1), int(at2), Chem.BondType.TRIPLE)
        elif bt == 5:
            new_mol.AddBond(int(at1), int(at2), Chem.BondType.SINGLE)
        elif bt == 6:
            new_mol.AddBond(int(at1), int(at2), Chem.BondType.DOUBLE)
        else:
            print(f"warning {bt=}")
            new_mol.AddBond(int(at1), int(at2), Chem.BondType.SINGLE)

    new_mol = new_mol.GetMol()
    # Chem.SanitizeMol(new_mol)
    new_mol = _add_conformer_to_mol(
        new_mol, np.array(prot_array.coord, dtype=np.float64)
    )
    # 去掉原始cif提供的H，感觉不太靠谱
    new_mol = Chem.RemoveAllHs(new_mol)
    return new_mol


def _lig_atomarray_to_mol(lig_array: AtomArray, lig_smiles: str) -> Chem.Mol:
    # new_mol = Chem.RWMol()
    mol = Chem.MolFromSmiles(lig_smiles)
    new_mol = Chem.RemoveAllHs(mol)
    assert len(lig_array.coord) == len(new_mol.GetAtoms())

    new_mol = _add_conformer_to_mol(
        new_mol, np.array(lig_array.coord, dtype=np.float64)
    )
    new_mol = Chem.AddHs(new_mol, addCoords=True)
    return new_mol


def _get_table(block: str, key: str) -> pd.DataFrame:
    lines = block.split("\n")
    headers = []
    table = []
    for line in lines[1:-2]:
        line = line.strip()
        if line.startswith(f"_{key}"):
            headers.append(line.split(".")[-1])
        else:
            values = line.split()
            row = dict(zip(headers, values))
            table.append(row)
    df = pd.DataFrame(table)
    return df


def get_SEQRES_block(df_seq: pd.DataFrame, modified_res: dict = None) -> Dict[str, str]:
    seqres_block = dict()
    for chain_id, df_group in df_seq.groupby(by="pdb_strand_id"):
        seq_lines = []
        grp_num = len(df_group.index)
        block_cnt = 0
        row_AAs = []
        for i, (ridx, row) in enumerate(df_group.iterrows()):
            if i % 13 == 0:
                block_cnt += 1
                row_AAs = []

            if modified_res and row["mon_id"] in modified_res:
                row["mon_id"] = modified_res[row["mon_id"]]

            row_AAs.append(row["mon_id"])
            if (i + 1) % 13 == 0:
                seq_prefix = f"SEQRES{block_cnt:>4d}{chain_id:>2s}{grp_num:>5d}  "
                seq_prefix += " ".join(row_AAs)
                seq_lines.append(seq_prefix + "\n")

        if len(row_AAs) > 0:
            seq_prefix = f"SEQRES{block_cnt:>4d}{chain_id:>2s}{grp_num:>5d}  "
            seq_prefix += " ".join(row_AAs)
            seq_lines.append(seq_prefix + "\n")
            row_AAs = []

        seqres_block[chain_id] = "".join(seq_lines)

    # for key, value in seqres_block.items():
    #     print(value)

    return seqres_block


def get_modified_res(df_mod: pd.DataFrame) -> Dict[str, str]:
    modified_res = dict()
    for ridx, row in df_mod.iterrows():
        label_comp_id = row["label_comp_id"]
        modified_residue_id = row["modified_residue_id"]
        # modified_res_name = row['modified_residue_name']
        modified_res[label_comp_id] = modified_residue_id

    return modified_res


def extract_cif_metainfo(cif_reader: pdbx.CIFFile) -> None:
    df_seq = None
    df_modification = None
    if "pdbx_poly_seq_scheme" in cif_reader.block._categories:
        df_seq = _get_table(
            cif_reader.block._categories["pdbx_poly_seq_scheme"], "pdbx_poly_seq_scheme"
        )

    if "pdbx_modification_feature" in cif_reader.block._categories:
        df_modification = _get_table(
            cif_reader.block._categories["pdbx_modification_feature"],
            "pdbx_modification_feature",
        )

    seq_blocks, modified_res = None, None

    if df_modification is not None:
        modified_res = get_modified_res(df_modification)

    if df_seq is not None:
        seq_blocks = get_SEQRES_block(df_seq, modified_res)

    return seq_blocks, modified_res


def update_het_atom_array(
    this_array: AtomArray, modified_res: Dict[str, str]
) -> AtomArray:
    for res1, res2 in modified_res.items():
        this_array.hetero[this_array.res_name == res1] = False
    return this_array


def process_one_cif(
    model_cif_fn: Path, input_json_dir: Path, out_dir: Path, platform: str = "CPU:16"
):
    cif_reader = pdbx.CIFFile.read(str(model_cif_fn))
    this_array = pdbx.get_structure(
        cif_reader,
        model=1,
        use_author_fields=True,
        include_bonds=True,
        extra_fields=["charge"],
    )

    item_name = model_cif_fn.stem
    pdb_id, lig_ccd = [x.upper() for x in item_name.split("_")[:2]]
    input_json_fn = input_json_dir / f"{pdb_id}_{lig_ccd}.json"
    with open(input_json_fn, "r") as f:
        input_json = json.load(f)

    lig_smiles = None
    lig_name = None
    for seq_info in input_json["sequences"]:
        if "ligand" in seq_info:
            lig_smiles = seq_info["ligand"]["smiles"]
            lig_name = "LIG_" + seq_info["ligand"]["id"]
            break
    assert lig_smiles is not None
    assert lig_name is not None

    prot_array, lig_array, lig_chain_id = extract_protein_and_ligand_array(
        this_array, lig_name
    )

    prot_mol = _prot_atomarray_to_mol(prot_array)
    lig_mol = _lig_atomarray_to_mol(lig_array, lig_smiles)

    lig_mol.SetProp("_Name", f"{item_name}_{lig_chain_id}")

    cur_dir = out_dir / item_name
    cur_dir.mkdir(parents=True, exist_ok=True)

    # out_mae_fn = cur_dir / f"{item_name}_protein_addh.mae"
    out_pdb_fn = cur_dir / f"{item_name}_protein_nofix.pdb"
    out_sdf_fn = cur_dir / f"{item_name}_ligand_fix.sdf"

    # mol_to_file(prot_mol, out_mae_fn)
    rdmol_to_file(prot_mol, out_pdb_fn)
    rdmol_to_file(lig_mol, out_sdf_fn)

    if out_pdb_fn.exists():
        protein_mol, missing_residues = fix_pdb(out_pdb_fn)
        out_pdb_fn = cur_dir / f"{item_name}_protein_fix.pdb"
        out_mae_fn = cur_dir / f"{item_name}_protein_fix.mae"
        rdmol_to_file(protein_mol, out_pdb_fn)
        rdmol_to_file(protein_mol, out_mae_fn)

        relax_tool = ProLigRelax(
            protein_mol,
            missing_residues=missing_residues,
            platform=platform,
            ligand_ff="openff",
            is_constrain=False,
            # charge_name= 'am1bcc'
        )
        prot_mol_relax, lig_mol_relax = relax_tool.prepare_one_cplx_and_relax(
            out_sdf_fn
        )

        logger.info(f"Relax success: {model_cif_fn}")
        out_mae_fn_relax = cur_dir / f"{item_name}_protein_fix_relax.mae"
        rdmol_to_file(prot_mol_relax, out_mae_fn_relax)

        out_sdf_fn_relax = cur_dir / f"{item_name}_ligand_fix_relax.sdf"
        rdmol_to_file(lig_mol_relax, out_sdf_fn_relax)

        from openmm.app import PDBFile

        out_pdb_fn_relax = cur_dir / f"{item_name}_protein_fix_relax.pdb"
        with open(out_pdb_fn_relax, "w") as f:
            PDBFile.writeFile(
                relax_tool.receptor_omm.topology,
                relax_tool.receptor_omm.positions,
                file=f,
                keepIds=True,
            )

    return


@click.command()
@click.option("-idir", "--cif_dir", type=str, default="posex/alphafold3_data/af3_outputs")
@click.option(
    "-jdir", "--json_input_dir", type=str, default="posex/alphafold3_data/input"
)
@click.option(
    "-odir", "--out_dir", type=str, default="posex/alphafold3_data/af3_outputs_relaxed"
)
def run(**kwargs):
    cif_dir = Path(kwargs["cif_dir"])
    json_input_dir = Path(kwargs["json_input_dir"])
    out_dir = Path(kwargs["out_dir"])
    cif_fns = list(cif_dir.glob("*.cif"))
    for cif_fn in tqdm(cif_fns, desc='relaxing', total = len(cif_fns)):
        if (out_dir / f"{cif_fn.stem}/{cif_fn.stem}_protein_fix_relax.pdb").exists():
            continue
        logger.info(f'{"="*20} {cif_fn} {"="*20}')
        try:
            process_one_cif(cif_fn, json_input_dir, out_dir, platform="CUDA:0")

        except:
            tb = traceback.format_exc()
            with open('error.log', 'a') as f:
                f.write(f'{"="}*20\n{cif_fn}\n{tb}\n')
                
        # if cif_fn.stem == "7qta_uri_model":
        #     process_one_cif(cif_fn, json_input_dir, out_dir, platform='CUDA:0')


if __name__ == "__main__":
    run()
    # fix_cif('/hpc-cache-pfs/home/xyj/code/protein_ligand_docking_benchmark/posex/mmcif_raw/8UCB_X1T.cif')
