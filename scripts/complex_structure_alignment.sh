#! /bin/bash

# Check if both arguments are provided
if [ $# -ne 3 ]; then
    echo "Error: Missing arguments"
    echo "Usage: $0 <dataset> <model_type> <relax_mode>"
    echo "Example: $0 posex_self_dock alphafold3 false"
    exit 1
fi

# Get the dataset and model_type from command-line arguments
DATASET=$1
MODEL_TYPE=$2
RELAX_MODE=$3

# Set dataset folder based on DATASET
if [ "$DATASET" = "posex_self_dock" ]; then
    DATASET_FOLDER="data/dataset/posex/posex_self_docking_set"
elif [ "$DATASET" = "posex_cross_dock" ]; then
    DATASET_FOLDER="data/dataset/posex/posex_cross_docking_set"
elif [ "$DATASET" = "posex_supp" ]; then
    DATASET_FOLDER="data/dataset/posex/posex_supp_set"
elif [ "$DATASET" = "astex" ]; then
    DATASET_FOLDER="data/dataset/posex/astex_diverse_set"
else
    echo "Error: Unknown dataset ${DATASET}"
    exit 1
fi

if [ "$RELAX_MODE" = "true" ]; then
    MODEL_OUTPUT_FOLDER="data/benchmark/${DATASET}/${MODEL_TYPE}/processed"
elif [ "$RELAX_MODE" = "false" ]; then
    MODEL_OUTPUT_FOLDER="data/benchmark/${DATASET}/${MODEL_TYPE}/output"
else
    echo "Error: Unknown relax_mode ${RELAX_MODE}"
    exit 1
fi

python scripts/complex_structure_alignment.py \
    --input_file data/benchmark/${DATASET}/${DATASET}_benchmark.csv \
    --dataset_folder ${DATASET_FOLDER} \
    --model_output_folder ${MODEL_OUTPUT_FOLDER} \
    --model_type ${MODEL_TYPE}
