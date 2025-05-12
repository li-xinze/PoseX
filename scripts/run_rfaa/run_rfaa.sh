# !/bin/bash


# Check if model type is provided as argument
if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"

GPU_ID=2
RFAA_INPUT_FOLDER="data/benchmark/${DATASET}/rfaa/input"
RFAA_OUTPUT_FOLDER="data/benchmark/${DATASET}/rfaa/output"

# Create output folder
mkdir -p ${RFAA_OUTPUT_FOLDER}

# init conda
eval "$(conda shell.bash hook)"
conda activate rfaa

RFAA_REPO_FOLDER="/home/jiangyize/workspace/github/RoseTTAFold-All-Atom/"

for yaml_file in ${RFAA_INPUT_FOLDER}/*.yaml; do
    filename=$(basename "${yaml_file}" .yaml)
    cp ${yaml_file} ${RFAA_REPO_FOLDER}/rf2aa/config/inference/

    pushd ${RFAA_REPO_FOLDER}
    echo "Predicting ${filename} ..."
    CUDA_VISIBLE_DEVICES=${GPU_ID} python -m rf2aa.run_inference --config-name ${filename}
    popd

    echo "Copying predictions of ${filename} ..."
    mv ${RFAA_REPO_FOLDER}/predictions/${filename} ${RFAA_OUTPUT_FOLDER}/
    rm ${RFAA_REPO_FOLDER}/rf2aa/config/inference/${filename}.yaml
done

conda activate posex