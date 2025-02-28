#! /bin/bash

if [ $# -eq 0 ]; then
    echo "Error: Please provide dataset as argument"
    echo "Usage: $0 <dataset>"
    exit 1
fi

DATASET="$1"
NEURALPLEXER_INPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/neuralplexer/input"
NEURALPLEXER_OUTPUT_FOLDER="${PWD}/data/benchmark/${DATASET}/neuralplexer/output"
RUNNING_SCRIPTS="${PWD}/scripts/run_neuralplexer/run_neuralplexer.py"
MODEL_CHECKPOINT="${PWD}/data/microcyto/astex/neuralplexer/input/complex_structure_prediction.ckpt"

if [ ! -e $MODEL_CHECKPOINT ]; then
  echo "File path ${MODEL_CHECKPOINT} does not exist"
  echo "Please download the model checkpoint"
  exit 1
fi

docker run -it \
    --privileged=true \
    --gpus "device=0" \
    -v $NEURALPLEXER_INPUT_FOLDER:/input \
    -v $MODEL_CHECKPOINT:/input/complex_structure_prediction.ckpt \
    -v $NEURALPLEXER_OUTPUT_FOLDER:/output \
    -v $RUNNING_SCRIPTS:/run_neuralplexer.py \
    neuralplexer:latest \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate NeuralPLexer && python /run_neuralplexer.py \
    --input_folder /input \
    --output_folder /output \
    --model_checkpoint /input/complex_structure_prediction.ckpt"
