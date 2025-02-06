import os

if __name__ == "__main__":
    input_dir ="/home/hanj/protein_ligand_docking_benchmark/data/microcyto/posebusters/gnina/output/"
    for d in os.listdir(input_dir):
        ligand_path = os.path.join(input_dir, d, f"{d}_model_ligand.sdf")
        output_path = os.path.join(input_dir, d, f"{d}_ligand.sdf")
        os.system(f"cp {ligand_path} {output_path}")
        os.system(f"rm {ligand_path}")
