#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
TANKBIND_EXEC_FOLDER="path/to/TankBind"
TANKBIND_INPUT_FOLDER="data/microcyto/${DATASET}/tankbind/input"
TANKBIND_OUTPUT_FOLDER="data/microcyto/${DATASET}/tankbind/output"
GPU_ID=0

# init conda
eval "$(conda shell.bash hook)"
conda activate tankbind_py38

p2rank_exec_path="${TANKBIND_EXEC_FOLDER}/package/p2rank/prank"
python scripts/run_tankbind/run_tankbind.py \
    --tankbind_exec_dir ${TANKBIND_EXEC_FOLDER} \
    --p2rank_exec_path ${p2rank_exec_path} \
    --input_dir ${TANKBIND_INPUT_FOLDER} \
    --output_dir ${TANKBIND_OUTPUT_FOLDER} \
    --gpu_id ${GPU_ID}

conda activate microcyto