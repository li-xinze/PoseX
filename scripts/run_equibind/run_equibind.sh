#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
EQUIBIND_EXEC_FOLDER="path/to/equibind"
EQUIBIND_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/equibind/input"
EQUIBIND_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/equibind/output"

# init conda
eval "$(conda shell.bash hook)"
conda activate equibind

start_time=$(date +%s)
python scripts/run_equibind/run_equibind.py \
    --input_dir ${EQUIBIND_INPUT_FOLDER} \
    --output_dir ${EQUIBIND_OUTPUT_FOLDER} \
    --equibind_exec_dir ${EQUIBIND_EXEC_FOLDER} \
    --yml_path "configs_clean/inference.yml"
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Running time for ${DATASET}: ${cost_time}"