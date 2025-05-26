#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
DYNAMICBIND_EXEC_FOLDER="path/to/DynamicBind"
DYNAMICBIND_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/dynamicbind/input")
DYNAMICBIND_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/dynamicbind/output")
GPU_ID=1

# init conda
eval "$(conda shell.bash hook)"
conda activate dynamicbind

for protein_filepath in ${DYNAMICBIND_INPUT_FOLDER}/*.pdb; do
    itemname=$(basename "${protein_filepath}" .pdb)
    output_folder="${DYNAMICBIND_OUTPUT_FOLDER}/${itemname}"
    ligand_filepath="${DYNAMICBIND_INPUT_FOLDER}/${itemname}.csv"
    echo "Predicting ${itemname}..."
    python scripts/run_dynamicbind/run_dynamicbind.py \
        --itemname ${itemname} \
        --protein_filepath ${protein_filepath} \
        --ligand_filepath ${ligand_filepath} \
        --dynamicbind_exec_dir ${DYNAMICBIND_EXEC_FOLDER} \
        --output_dir ${output_folder} \
        --gpu_id ${GPU_ID}
done
