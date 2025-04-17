#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
INTERFORMER_EXEC_FOLDER="path/to/interformer"
INTERFORMER_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/interformer/input"
INTERFORMER_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/interformer/output"

# init conda
eval "$(conda shell.bash hook)"
conda activate interformer

export CUDA_VISIBLE_DEVICES="1"

start_time=$(date +%s)
python scripts/run_interformer/run_interformer.py \
    --input_dir ${INTERFORMER_INPUT_FOLDER} \
    --output_dir ${INTERFORMER_OUTPUT_FOLDER} \
    --interformer_exec_dir ${INTERFORMER_EXEC_FOLDER}
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Running time for ${DATASET}: ${cost_time}"