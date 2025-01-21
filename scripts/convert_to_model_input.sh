#! /bin/bash


# Check if both arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Missing arguments"
    echo "Usage: $0 <dataset> <model_type>"
    echo "Example: $0 astex alphafold3"
    exit 1
fi

# Get the dataset and model_type from command-line arguments
DATASET=$1
MODEL_TYPE=$2

python scripts/convert_to_model_input.py \
    --input_file data/benchmark/${DATASET}/${DATASET}_benchmark.csv \
    --output_folder data/benchmark/${DATASET}/${MODEL_TYPE}/input \
    --model_type ${MODEL_TYPE}
