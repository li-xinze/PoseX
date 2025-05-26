#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
FABIND_EXEC_FOLDER="path/to/fabind"
CKPT_PATH="${FABIND_EXEC_FOLDER}/ckpt/best_model.bin"
FABIND_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/fabind/input")
FABIND_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/fabind/output")
GPU_ID=0

# init conda
eval "$(conda shell.bash hook)"
conda activate fabind

input_csv_path="${FABIND_INPUT_FOLDER}/ligand.csv"
input_data_dir="${FABIND_INPUT_FOLDER}/protein"

python scripts/run_fabind/run_fabind.py \
    --input_csv_path ${input_csv_path} \
    --input_data_dir ${input_data_dir} \
    --output_dir ${FABIND_OUTPUT_FOLDER} \
    --fabind_exec_dir ${FABIND_EXEC_FOLDER} \
    --ckpt_path ${CKPT_PATH} \
    --gpu_id ${GPU_ID}
