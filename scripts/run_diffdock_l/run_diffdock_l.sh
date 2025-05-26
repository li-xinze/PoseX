#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
DIFFDOCK_EXEC_FOLDER="path/to/DiffDock_L"
DIFFDOCK_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/diffdock_l/input")
DIFFDOCK_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/diffdock_l/output")


# init conda
eval "$(conda shell.bash hook)"
conda activate diffdock

input_csv_path="${DIFFDOCK_INPUT_FOLDER}/data.csv"
config_path="${DIFFDOCK_EXEC_FOLDER}/default_inference_args.yaml"
python scripts/run_diffdock_l/run_diffdock_l.py \
    --input_csv_path ${input_csv_path} \
    --config_path ${config_path} \
    --output_dir ${DIFFDOCK_OUTPUT_FOLDER} \
    --diffdock_exec_dir ${DIFFDOCK_EXEC_FOLDER}
