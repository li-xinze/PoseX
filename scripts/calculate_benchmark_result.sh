#! /bin/bash

# Check if both arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Missing arguments"
    echo "Usage: $0 <dataset> <model_type>"
    echo "Example: $0 posex_self_dock alphafold3"
    exit 1
fi

# Get the dataset and model_type from command-line arguments
DATASET=$1
MODEL_TYPE=$2

# Set dataset folder based on DATASET
if [ "$DATASET" = "posex_self_dock" ]; then
    DATASET_FOLDER="data/dataset/posex/posex_self_docking_set"
elif [ "$DATASET" = "posex_cross_dock" ]; then
    DATASET_FOLDER="data/dataset/posex/posex_cross_docking_set"
else
    echo "Error: Unknown dataset ${DATASET}"
    exit 1
fi


python scripts/calculate_benchmark_result.py \
    --input_file data/benchmark/${DATASET}/${DATASET}_benchmark.csv \
    --dataset_folder ${DATASET_FOLDER} \
    --model_output_folder data/benchmark/${DATASET}/${MODEL_TYPE}/output \
    --model_type ${MODEL_TYPE} \
    --dataset ${DATASET}
