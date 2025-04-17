#! /bin/bash

if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
DEEPDOCK_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/deepdock/input"
DEEPDOCK_RUNNING_FOLDER="${PWD}/scripts/run_deepdock"

EVALUATE_FILE="scripts/run_deepdock/evaluate.py"
for subdir in $DEEPDOCK_INPUT_FOLDER/*; do
    output_dir=$subdir
    cp $EVALUATE_FILE $output_dir
done

# DeepDock only accept mol2 input for ligand
# Install openbabel to convert sdf file into mol2
echo "convert sdf file to mol2 using obabel"
eval "$(conda shell.bash hook)"
conda activate openbabel

for subdir in $DEEPDOCK_INPUT_FOLDER/*; do
    subdir_name=$(basename $subdir)
    sdf_file="$subdir/"$subdir_name"_ligand.sdf"
    output_file="$subdir/"$subdir_name"_ligand.mol2"
    obabel $sdf_file -O $output_file
done

start_time=$(date +%s)
docker run -it \
    -v $DEEPDOCK_INPUT_FOLDER:/DeepDock/eval \
    -v $DEEPDOCK_RUNNING_FOLDER:/DeepDock/run \
    omendezlucio/deepdock \
    python DeepDock/run/run_deepdock.py --input_folder /DeepDock/eval
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "Running time for ${DATASET}: ${cost_time}"