import multiprocessing as mp
import traceback
from functools import partial
from pathlib import Path
from typing import *
import click
from tqdm import tqdm
import os

from dataset.utils.pdb_helper import load_amber_xml
from dataset.utils.pdb_process import (
    fixer_into_protein_mol,
    pdb_to_fixer,
    protein_mol_to_file,
    ligand_rdmol_to_file,
)

from dataset.utils.openmm_helper import ProLigRelax
from dataset.utils.common_helper import create_logger


logger = create_logger(__name__)
NUM_CPUS = mp.cpu_count()


def run_once(
    pdb_fn: Path,
    output_dir: Path,
    residues_tables: tuple[Dict, Dict, Dict],
    num_threads: int = 1,
):
    print(pdb_fn)
    tmp_items = pdb_fn.stem.split("_")
    tmp_name = "_".join(tmp_items[:2])
    sdf_fn = pdb_fn.parent / f"{tmp_name}_model_ligand.sdf"
    cif_fn = pdb_fn.parent / f"{tmp_items[0]}.cif"
    work_dir = output_dir / f"{tmp_name}"
    work_dir.mkdir(parents=True, exist_ok=True)
    out_pdb_fn = work_dir / f"{tmp_name}_protein_step1.pdb"

    out_pdb_fn2 = work_dir / f"{tmp_name}_protein_step2.pdb"
    out_ligand_fn2 = work_dir / f"{tmp_name}_ligand_step2.sdf"
    if out_ligand_fn2.exists() and out_pdb_fn2.exists():
        return True

    processed_dir = output_dir / pdb_fn.parent.name
    processed_dir.mkdir(parents=True, exist_ok=True)
    try:

        fixer_noh = pdb_to_fixer(pdb_fn, cif_fn)
        prot_mol_h = fixer_into_protein_mol(fixer_noh, residues_tables)
        protein_mol_to_file(prot_mol_h, out_pdb_fn)
        print(f"{out_pdb_fn=}")
        # relax
        logger.info(f"Run relaxing...")
        relax_tool = ProLigRelax(
            prot_mol_h,
            missing_residues=[],
            platform=f"CPU:{num_threads}",
            ligand_ff="openff",
            charge_name="mmff94",
            is_restrain=(True, "main"),
        )
        receptor_rdmol_relaxed, ligand_rdmol_relaxed = (
            relax_tool.prepare_one_cplx_and_relax(sdf_fn)
        )

        logger.info(f"Relax success: {cif_fn=}")

        ligand_rdmol_to_file(ligand_rdmol_relaxed, out_ligand_fn2)
        logger.info(f"Writing ligand to {out_ligand_fn2}")

        protein_mol_to_file(receptor_rdmol_relaxed, out_pdb_fn2)
        logger.info(f"Writing protein to {out_pdb_fn2}")
        logger.info(f'{"#"*10} Succeed: {pdb_fn=}')

        del relax_tool

    except:
        tb = traceback.format_exc()
        tmp_error_log = output_dir.parent / "error.log"
        with open(tmp_error_log, "a") as f:
            f.write(f"#PDBNAME#: {tmp_name}\n{tb}")
        logger.info(f"{tmp_error_log=}")

    return True


def run_batch(pdb_fns, output_dir, num_proc: int = 6):
    residues_tables = load_amber_xml()
    num_threads = int(NUM_CPUS / num_proc)
    with mp.Pool(num_proc) as pool:
        func = partial(
            run_once,
            output_dir=output_dir,
            residues_tables=residues_tables,
            num_threads=num_threads,
        )
        results = list(
            tqdm(
                pool.imap(func, pdb_fns, chunksize=1),
                total=len(pdb_fns),
                desc="Processing",
            )
        )

    print(f"{sum(results)=}")


@click.command()
@click.option(
    "--input_dir",
    type=str,
    default="{path_to_data}/fabind/posex_self_dock/output",
)
@click.option("--output_dir", type=str, default=None)
@click.option("--num_proc", type=int, default=1)
def main(**kwargs):
    input_dir = Path(kwargs["input_dir"])
    if kwargs["output_dir"] is None:
        output_dir = input_dir.parent / "processed"
        error_fn = input_dir.parent / "error.log"
        if error_fn.exists():
            error_fn.unlink(missing_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path(kwargs["output_dir"])
    pdb_fns = list(input_dir.glob("*/*.pdb"))

    run_batch(pdb_fns, output_dir, int(kwargs["num_proc"]))


if __name__ == "__main__":
    main()
