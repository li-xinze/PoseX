#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
UNIMOL_EXEC_FOLDER="/home/hanj/Docking/Uni-Mol/unimol_docking_v2/"
CKPT_PATH="${UNIMOL_EXEC_FOLDER}/ckpt/unimol_docking_v2_240517.pt"
UNIMOL_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/unimol/input"
UNIMOL_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/unimol/output"

# init conda
eval "$(conda shell.bash hook)"
conda activate unicore

python scripts/run_unimol/run_unimol.py \
    --input_dir ${UNIMOL_INPUT_FOLDER} \
    --output_dir ${UNIMOL_OUTPUT_FOLDER} \
    --unimol_exec_dir ${UNIMOL_EXEC_FOLDER} \
    --ckpt_path ${CKPT_PATH} \

#conda activate microcyto