#!/bin/bash

PART_ID=3

# Check if model type is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"

# Run AlphaFold3
docker run -it --rm --gpus all --shm-size=32g -e CUDA_VISIBLE_DEVICES=${PART_ID} \
	-v ./data/microcyto/${DATASET}/alphafold3/input/part_${PART_ID}:/root/af_input \
	-v ./data/microcyto/${DATASET}/alphafold3/output:/root/af_output \
	-v /data/dataset/alphafold3/models:/root/models \
	-v /data/dataset/alphafold3/databases:/root/public_databases \
	brandonsoubasis/alphafold3 \
	python run_alphafold.py \
    --input_dir=/root/af_input \
	--model_dir=/root/models \
	--output_dir=/root/af_output