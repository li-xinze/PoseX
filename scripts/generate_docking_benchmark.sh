#! /bin/bash


# Check if model type is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"

python scripts/generate_docking_benchmark.py \
    --input_folder data/posebusters \
    --output_folder data/microcyto/${DATASET} \
    --dataset ${DATASET}
