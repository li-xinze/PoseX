#! /bin/bash

BOLTZ_INPUT_FOLDER="data/microcyto/posebusters/blotz/input"
BOLTZ_OUTPUT_FOLDER="data/microcyto/posebusters/blotz/output"
GPU_ID=3

# init conda
eval "$(conda shell.bash hook)"
conda activate boltz

for yaml_file in ${BOLTZ_INPUT_FOLDER}/*.yaml; do
    filename=$(basename "${yaml_file}" .yaml)
    output_folder="${BOLTZ_OUTPUT_FOLDER}/${filename}"

    echo "Predicting ${yaml_file} to ${output_folder} ..."
    CUDA_VISIBLE_DEVICES=${GPU_ID} boltz predict ${yaml_file} --out_dir ${output_folder} --cache /data/models/boltz --use_msa_server
done

conda activate microcyto