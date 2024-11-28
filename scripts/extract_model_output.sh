#! /bin/bash


# MODEL_TYPE="alphafold3"
MODEL_TYPE="chai"

python scripts/extract_model_output.py \
    --input_file data/microcyto/posebusters/posebusters_benchmark.csv \
    --output_folder data/microcyto/posebusters/${MODEL_TYPE}/output \
    --model_type ${MODEL_TYPE}
