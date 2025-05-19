
<div align="center">
  <img src="figures/logo.png" width="400"/>
</div>
<div>&nbsp;</div>

# PoseX: AI Defeats Physics-based Methods on Protein Ligand Cross-Docking

<div align="center">

[![arXiv](https://img.shields.io/badge/arXiv-2505.01700-b31b1b.svg)](https://arxiv.org/abs/2505.01700)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[![Dataset on HF](https://huggingface.co/datasets/huggingface/badges/resolve/main/dataset-on-hf-sm.svg)](https://huggingface.co/datasets/CataAI/PoseX)
</div>

PoseX is a comprehensive benchmark dataset designed to evaluate molecular docking algorithms for predicting protein-ligand binding poses. It includes the construction process of Self-Docking and Cross-Docking datasets, as well as complete evaluation codes for different docking tools.

[Online Leaderboard](http://116.62.156.219)

[Dataset Repository](https://huggingface.co/datasets/CataAI/PoseX)

---
<div align="center" style="font-size:24px;">
PoseX Self-Docking Result
</div>
<div align="center">
  <img src="figures/posex_self_dock.png"/>
</div>
<div>&nbsp;</div>

<div align="center" style="font-size:24px;">
PoseX Cross-Docking Result
</div>
<div align="center">
  <img src="figures/posex_cross_dock.png"/>
</div>
<div>&nbsp;</div>

## Contents
- [PoseX: AI Defeats Physics-based Methods on Protein Ligand Cross-Docking](#posex-ai-defeats-physics-based-methods-on-protein-ligand-cross-docking)
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
  - [Citations](#citations)

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



## Citations

If you are interested in our work or use our data and code, please cite the following article:

```bibtex
@misc{jiang2025posexaidefeatsphysics,
      title={PoseX: AI Defeats Physics Approaches on Protein-Ligand Cross Docking}, 
      author={Yize Jiang and Xinze Li and Yuanyuan Zhang and Jin Han and Youjun Xu and Ayush Pandit and Zaixi Zhang and Mengdi Wang and Mengyang Wang and Chong Liu and Guang Yang and Yejin Choi and Wu-Jun Li and Tianfan Fu and Fang Wu and Junhong Liu},
      year={2025},
      eprint={2505.01700},
      archivePrefix={arXiv},
      primaryClass={cs.LG},
      url={https://arxiv.org/abs/2505.01700}, 
}
```