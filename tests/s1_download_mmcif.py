from pathlib import Path
import wget

in_dir = "posex/posex_self_docking_set"
out_dir = "posex/mmcif_raw"


def download_one(pdb_id: str, save_fn: Path):
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"  # 请替换为实际的URL
    # 下载文件
    wget.download(url, str(save_fn))


def run():
    input_fns = list(Path(in_dir).glob("*/*.json"))
    cur_dir = Path(out_dir)
    cur_dir.mkdir(exist_ok=True, parents=True)
    for input_fn in input_fns[:10]:
        item_name = input_fn.parent.stem
        pdb_id = item_name.split("_")[0].lower()
        save_fn = cur_dir / f"{item_name}.cif"
        download_one(pdb_id, Path(save_fn))
        print(input_fn)


run()
