#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
#INTERFORMER_EXEC_FOLDER="path/to/interformer"
INTERFORMER_EXEC_FOLDER="/home/hanj/Docking/Interformer/"
INTERFORMER_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/interformer/input"
INTERFORMER_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/interformer/output"

# init conda
eval "$(conda shell.bash hook)"
conda activate interformer

python scripts/run_interformer/run_interformer.py \
    --input_dir ${INTERFORMER_INPUT_FOLDER} \
    --output_dir ${INTERFORMER_OUTPUT_FOLDER} \
    --interformer_exec_dir ${INTERFORMER_EXEC_FOLDER} \

#conda activate microcyto