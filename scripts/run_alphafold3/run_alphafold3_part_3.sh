#!/bin/bash

PART_ID=3

# Run AlphaFold3
docker run -it --rm --gpus all --shm-size=32g -e CUDA_VISIBLE_DEVICES=${PART_ID} \
	-v /home/jiangyize/workspace/projects/protein_ligand_docking_benchmark/data/microcyto/posebusters/alphafold3/input/part_${PART_ID}:/root/af_input \
	-v /home/jiangyize/workspace/projects/protein_ligand_docking_benchmark/data/microcyto/posebusters/alphafold3/output:/root/af_output \
	-v /data/dataset/alphafold3/models:/root/models \
	-v /data/dataset/alphafold3/databases:/root/public_databases \
	brandonsoubasis/alphafold3 \
	python run_alphafold.py \
    --input_dir=/root/af_input \
	--model_dir=/root/models \
	--output_dir=/root/af_output
