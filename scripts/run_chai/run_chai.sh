#! /bin/bash


# Check if model type is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"

CHAI_INPUT_FOLDER="data/microcyto/${DATASET}/chai/input"
CHAI_OUTPUT_FOLDER="data/microcyto/${DATASET}/chai/output"
GPU_ID=2

# init conda
eval "$(conda shell.bash hook)"
conda activate chai

for fasta_file in ${CHAI_INPUT_FOLDER}/*.fasta; do
    filename=$(basename "${fasta_file}" .fasta)
    output_folder="${CHAI_OUTPUT_FOLDER}/${filename}"

    echo "Predicting ${fasta_file} to ${output_folder} ..."
    python scripts/run_chai/run_chai.py \
        --fasta_file ${fasta_file} \
        --output_dir ${output_folder} \
        --gpu_id ${GPU_ID}
done

conda activate microcyto
