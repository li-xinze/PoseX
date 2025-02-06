#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
#INTERFORMER_EXEC_FOLDER="path/to/equibind"
EQUIBIND_EXEC_FOLDER="/home/hanj/Docking/EquiBind/"
EQUIBIND_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/equibind/input"
EQUIBIND_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/equibind/output"

# init conda
eval "$(conda shell.bash hook)"
conda activate equibind

python scripts/run_equibind/run_equibind.py \
    --input_dir ${EQUIBIND_INPUT_FOLDER} \
    --output_dir ${EQUIBIND_OUTPUT_FOLDER} \
    --equibind_exec_dir ${EQUIBIND_EXEC_FOLDER} \
    --yml_path "configs_clean/inference.yml"

#conda activate microcyto