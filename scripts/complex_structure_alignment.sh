#! /bin/bash

# Check if both arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Missing arguments"
    echo "Usage: $0 <dataset> <model_type>"
    echo "Example: $0 posebusters alphafold3"
    exit 1
fi

# Get the dataset and model_type from command-line arguments
DATASET=$1
MODEL_TYPE=$2

# Set dataset folder based on DATASET
if [ "$DATASET" = "posebusters" ]; then
    DATASET_FOLDER="data/posebusters/posebusters_benchmark_set"
elif [ "$DATASET" = "astex" ]; then
    DATASET_FOLDER="data/posebusters/astex_diverse_set"
else
    echo "Error: Unknown dataset ${DATASET}"
    exit 1
fi

python scripts/complex_structure_alignment.py \
    --input_file data/microcyto/${DATASET}/${DATASET}_benchmark.csv \
    --dataset_folder ${DATASET_FOLDER} \
    --model_output_folder data/microcyto/${DATASET}/${MODEL_TYPE}/output \
    --model_type ${MODEL_TYPE}
