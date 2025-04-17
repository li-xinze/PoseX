#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
UNIMOL_EXEC_FOLDER="path/to/unimol"
CKPT_PATH="${UNIMOL_EXEC_FOLDER}/ckpt/unimol_docking_v2_240517.pt"
UNIMOL_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/unimol/input"
UNIMOL_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/unimol/output"

# init conda
eval "$(conda shell.bash hook)"
conda activate unicore

export MKL_SERVICE_FORCE_INTEL=1
export CUDA_VISIBLE_DEVICES=1


start_time=$(date +%s)
python scripts/run_unimol/run_unimol.py \
    --input_dir ${UNIMOL_INPUT_FOLDER} \
    --output_dir ${UNIMOL_OUTPUT_FOLDER} \
    --unimol_exec_dir ${UNIMOL_EXEC_FOLDER} \
    --ckpt_path ${CKPT_PATH}
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Running time for ${DATASET}: ${cost_time}"