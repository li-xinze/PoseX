#! /bin/bash


# Check if dataset is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
DIFFDOCK_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/diffdock_pocket/input"
DIFFDOCK_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/diffdock_pocket/output"
DIFFDOCK_EXEC_PATH="path/to/Diffdock_pocket"

cd $DIFFDOCK_EXEC_PATH

start_time=$(date +%s)
python inference.py --protein_ligand_csv "${DIFFDOCK_INPUT_FOLDER}/example.csv" --out_dir "${DIFFDOCK_OUTPUT_FOLDER}/results" --batch_size 12 --samples_per_complex 40 --keep_local_structures
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Running time for ${DATASET}: ${cost_time}"