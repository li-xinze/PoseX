
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
    - [1. Generate Benchmark CSV Data](#1-generate-benchmark-csv-data)
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
cd PoseX
```
We recommend using `mamba` to manage the Python environment. For more information on how to install mamba, see [Miniforge](https://github.com/conda-forge/miniforge).
Once `mamba` is installed, we can run the following command to install the basic environment.
```bash
mamba create -f environments/base.yaml
mamba activate posex
```
For a specific molecular docking tool, we can use the corresponding environment file in the `environments` folder. Take `Chai-1` as an example:

```bash
pip install -r environments/chai-1.txt
```

## Benchmark Data
For information about creating a dataset from scratch, please refer to the data construction [README](./dataset/README.md).

## Benchmark Pipeline

This project provides a complete pipeline for:
- Generate the csv file required for the benchmark based on the `PoseX` dataset
- Converting docking data into model-specific input formats
- Running different docking models or tools
- Energy minimization of molecular docking results
- Extracting and aligning model outputs
- Calculating evaluation metrics using `PoseBusters`


### 1. Generate Benchmark CSV Data

Generate benchmark CSV files containing protein sequences, ligand SMILES, and other metadata:

```bash
bash ./scripts/generate_docking_benchmark.sh <dataset>

# --------------- Example --------------- #
# Astex
bash ./scripts/generate_docking_benchmark.sh astex
# PoseX Self-Docking
bash ./scripts/generate_docking_benchmark.sh posex_self_dock
# PoseX Cross-Docking
bash ./scripts/generate_docking_benchmark.sh posex_cross_dock
```

### 2. Convert to Model Inputs

Convert benchmark CSV files to model-specific input formats:

```bash
bash ./scripts/convert_to_model_input.sh <dataset> <model_type>

# --------------- Example --------------- #
# PoseX Self-Docking (AlphaFold3)
bash ./scripts/convert_to_model_input.sh posex_self_dock alphafold3
```

### 3. Run Docking Models

Run different docking models:

```bash
bash ./scripts/run_<model_type>/run_<model_type>.sh <dataset>

# --------------- Example --------------- #
# PoseX Self-Docking (Alphafold3)
bash ./scripts/run_alphafold3/run_alphafold3.sh posex_self_dock
```

### 4. Extract Model Outputs

Extract predicted structures from model outputs:

```bash
bash ./scripts/extract_model_output.sh <dataset> <model_type>

# --------------- Example --------------- #
# PoseX Self-Docking (AlphaFold3)
bash ./scripts/extract_model_output.sh posex_self_dock alphafold3
```

### 5. Energy Minimization

```bash
python -m scripts.relax_model_outputs --input_dir <input_dir> --cif_dir <cif_dir>
```

### 6. Align Predicted Structures

Align predicted structures to the reference structures:

```bash
bash ./scripts/complex_structure_alignment.sh <dataset> <model_type> <relax_mode>

# --------------- Example --------------- #
# PoseX Self-Docking (AlphaFold3) (Using Relax)
bash ./scripts/complex_structure_alignment.sh posex_self_dock alphafold3 true
```

### 7. Calculate Benchmark Result

Calculate evaluation metrics using PoseBusters:

```bash
bash ./scripts/calculate_benchmark_result.sh <dataset> <model_type> <relax_mode>

# --------------- Example --------------- #
# PoseX Self-Docking (AlphaFold3) (Using Relax)
bash ./scripts/calculate_benchmark_result.sh posex_self_dock alphafold3 true
```

## Acknowledgements

## License
- **Code**: Licensed under the [MIT License](https://opensource.org/licenses/MIT). 
- **Dataset**: Licensed under [Creative Commons Attribution 4.0 International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/). See [PoseX Dataset](https://huggingface.co/datasets/CataAI/PoseX) for details.



## Cite

If you are interested in our work or use our data and code, please cite the following article:
