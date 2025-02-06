#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
PROTENIX_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/protenix/input"
PROTENIX_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/protenix/output"
# init conda
eval "$(conda shell.bash hook)"
conda activate protenix

protenix predict --input $PROTENIX_INPUT_FOLDER --out_dir  $PROTENIX_OUTPUT_FOLDER --seeds 101 --use_msa_server

#conda activate microcyto