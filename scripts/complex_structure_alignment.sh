#! /bin/bash

# MODEL_TYPE="alphafold3"
MODEL_TYPE="chai"

python scripts/complex_structure_alignment.py \
    --input_file data/microcyto/posebusters/posebusters_benchmark.csv \
    --dataset_folder data/posebusters/posebusters_benchmark_set \
    --model_output_folder data/microcyto/posebusters/${MODEL_TYPE}/output \
    --model_type ${MODEL_TYPE}
