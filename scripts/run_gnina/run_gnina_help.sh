#!/bin/bash

INPUT_FOLDER="$1"
OUTPUT_FOLDER="$2"

for pdb_ccd in "$INPUT_FOLDER"/*; do
    if [ -d "$pdb_ccd" ]; then
        protein_path="$pdb_ccd/${pdb_ccd##*/}_protein.pdb"
        ligand_path="$pdb_ccd/${pdb_ccd##*/}_ligand_start_conf.sdf"

        output_dir="$OUTPUT_FOLDER/$(basename "$pdb_ccd")"
        output_path="$output_dir/$(basename "$pdb_ccd")_ligand.sdf"

        mkdir -p "$output_dir"

        cmd="gnina -r $protein_path -l $ligand_path --autobox_ligand $ligand_path -o $output_path"

        eval "$cmd"
    fi
done
