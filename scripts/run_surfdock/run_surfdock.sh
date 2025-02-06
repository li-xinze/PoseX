#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
SURFDOCK_INPUT_FOLDER="${PWD}/data/microcyto/${DATASET}/surfdock/input"
SURFDOCK_OUTPUT_FOLDER="${PWD}/data/microcyto/${DATASET}/surfdock/output"
RUNNING_SCRIPTS="${PWD}/scripts/run_surfdock/run_surfdock_help.sh"

docker run -it \
    --privileged=true \
    --gpus all \
    -e CUDA_VISIBLE_DEVICES="1" \
    -v /home/hanj/workplace/SurfDock/:/SurfDock \
    -v $SURFDOCK_INPUT_FOLDER:/input \
    -v $SURFDOCK_OUTPUT_FOLDER:/output \
    -v $RUNNING_SCRIPTS:/run_surfdock_help.sh \
    surfdock:v1 \
    bash -c "/run_surfdock_help.sh /input /output"
