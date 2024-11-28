#! /bin/bash

# MODEL_TYPE="alphafold3"
# MODEL_TYPE="chai"
MODEL_TYPE="blotz"

python scripts/convert_to_model_input.py \
    --input_file data/microcyto/posebusters/posebusters_benchmark.csv \
    --output_folder data/microcyto/posebusters/${MODEL_TYPE}/input \
    --model_type ${MODEL_TYPE}
