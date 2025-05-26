#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
DIFFDOCK_EXEC_FOLDER="path/to/DiffDock"
MODEL_DIR="${DIFFDOCK_EXEC_FOLDER}/workdir/paper_score_model"
CONFIDENCE_MODEL_DIR="${DIFFDOCK_EXEC_FOLDER}/workdir/paper_confidence_model"
DIFFDOCK_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/diffdock/input")
DIFFDOCK_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/diffdock/output")
GPU_ID=1

# init conda
eval "$(conda shell.bash hook)"
conda activate diffdock

input_csv_path="${DIFFDOCK_INPUT_FOLDER}/data.csv"
python scripts/run_diffdock/run_diffdock.py \
    --input_csv_path ${input_csv_path} \
    --output_dir ${DIFFDOCK_OUTPUT_FOLDER} \
    --gpu_id ${GPU_ID} \
    --diffdock_exec_dir ${DIFFDOCK_EXEC_FOLDER} \
    --model_dir ${MODEL_DIR} \
    --confidence_model_dir ${CONFIDENCE_MODEL_DIR}
