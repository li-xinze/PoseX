#! /bin/bash

if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
GNINA_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/gnina/input"
GNINA_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/gnina/output"
RUNNING_SCRIPTS="${PWD}/scripts/run_gnina/run_gnina_help.sh"

docker run -it \
    --privileged=true \
    --gpus all \
    -e CUDA_VISIBLE_DEVICES="1" \
    -v $GNINA_INPUT_FOLDER:/input \
    -v $GNINA_OUTPUT_FOLDER:/output \
    -v $RUNNING_SCRIPTS:/run_gnina_help.sh \
    gnina/gnina:latest \
    bash -c "/run_gnina_help.sh /input /output"

