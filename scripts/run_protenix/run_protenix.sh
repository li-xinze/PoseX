#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
PROTENIX_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/protenix/input"
PROTENIX_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/protenix/output"

eval "$(conda shell.bash hook)"
conda activate protenix

export CUDA_VISIBLE_DEVICES="1"

start_time=$(date +%s)
protenix predict --input $PROTENIX_INPUT_FOLDER --out_dir  $PROTENIX_OUTPUT_FOLDER --seeds 101 --use_msa_server
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Running time for ${DATASET}: ${cost_time}"

