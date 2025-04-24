
<div align="center">
  <img src="figures/logo.png" width="400"/>
</div>
<div>&nbsp;</div>

# PoseX: A Molecular Docking Benchmark

PoseX is a comprehensive benchmark dataset designed to evaluate molecular docking algorithms for predicting protein-ligand binding poses. It includes the construction process of Self-Docking and Cross-Docking datasets, as well as complete evaluation codes for different docking tools.

## Contents
- [PoseX: A Molecular Docking Benchmark](#posex-a-molecular-docking-benchmark)
  - [Contents](#contents)
  - [Installation](#installation)
  - [Benchmark Data](#benchmark-data)
  - [Benchmark Pipeline](#benchmark-pipeline)
    - [1. Generate Benchmark Dataset](#1-generate-benchmark-dataset)
    - [2. Convert to Model Inputs](#2-convert-to-model-inputs)
    - [3. Run Docking Models](#3-run-docking-models)
    - [4. Extract Model Outputs](#4-extract-model-outputs)
    - [5. Energy Minimization](#5-energy-minimization)
    - [6. Align Predicted Structures](#6-align-predicted-structures)
    - [7. Calculate Benchmark Result](#7-calculate-benchmark-result)
  - [Acknowledgements](#acknowledgements)
  - [License](#license)
  - [Cite](#cite)

## Installation

Install PoseX directly from Github to get the latest updates.
```bash
git clone https://github.com/CataAI/PoseX.git
```
We recommend using `mamba` to manage the Python environment. For more information on how to install mamba, see [Miniforge](https://github.com/conda-forge/miniforge).
Once `mamba` is installed, we can run the following command to install the basic environment.
```bash
mamba create -f environments/base.yaml
mamba activate posex
```
For a specific molecular docking tool, we can use the corresponding environment file in the `environments` folder. Take `Chai-1` as an example:

```bash
pip install -r environments/chai.txt
```

## Benchmark Data
For information about creating a dataset from scratch, please refer to the data construction [README](./dataset/README.md).

## Benchmark Pipeline

This project provides a complete pipeline for:
- Generating benchmark datasets from PoseX and Astex
- Converting data into model-specific input formats
- Running different docking models
- Energy minimization of molecular docking results
- Extracting and aligning model outputs
- Calculating evaluation metrics using PoseBusters


### 1. Generate Benchmark Dataset

Generate benchmark CSV files containing protein sequences, ligand SMILES, and other metadata:

```bash
bash ./scripts/generate_docking_benchmark.sh <dataset>
```

### 2. Convert to Model Inputs

Convert benchmark CSV files to model-specific input formats:

```bash
bash ./scripts/convert_to_model_input.sh <dataset> <model_type>
```

### 3. Run Docking Models

Run different docking models:

For example, run Boltz-1:
```bash
bash ./scripts/run_boltz/run_boltz.sh <dataset>
```

### 4. Extract Model Outputs

Extract predicted structures from model outputs:

```bash
bash ./scripts/extract_model_output.sh <dataset> <model_type>
```

### 5. Energy Minimization

```
python -m scripts.relax_model_outputs --input_dir <input_dir> --cif_dir <cif_dir>
```

### 6. Align Predicted Structures

Align predicted structures to the reference structures:

```bash
bash ./scripts/complex_structure_alignment.sh <dataset> <model_type> <relax_mode>
```

### 7. Calculate Benchmark Result

Calculate evaluation metrics using PoseBusters:

```bash
bash ./scripts/calculate_benchmark_result.sh <dataset> <model_type> <relax_mode>
```

## Acknowledgements

## License
- **Code**: Licensed under the [MIT License](https://opensource.org/licenses/MIT). 
- **Dataset**: Licensed under [Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/). See [PoseX Dataset](https://huggingface.co/datasets/CataAI/PoseX) for details.



## Cite

If you are interested in our work or use our data and code, please cite the following article:
