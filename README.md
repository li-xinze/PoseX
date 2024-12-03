# Docking Benchmark

A comprehensive benchmarking framework for evaluating protein-ligand docking models, including AlphaFold3, Chai-1, Boltz-1 etc.

## Project Overview

This project provides a complete pipeline for:
- Generating benchmark datasets from PoseBusters and Astex
- Converting data into model-specific input formats
- Running different docking models
- Extracting and aligning model outputs
- Calculating evaluation metrics using PoseBusters

## Directory Structure
```plaintext
.
├── data     # Data directory (gitignored)
│   ├── microcyto      # Processed benchmark data
│   │   ├── astex
│   │   └── posebusters
│   └── posebusters    # Raw PoseBusters dataset
├── environments    # Environment configuration files
└── scripts
    ├── calculate_benchmark_result.sh     # Calculate benchmark result
    ├── complex_structure_alignment.sh    # Align predicted structures
    ├── convert_to_model_input.sh         # Convert to model inputs
    ├── extract_model_output.sh           # Extract model outputs
    ├── generate_docking_benchmark.sh     # Generate benchmark datasets
    ├── run_alphafold3                    # Run AlphaFold3
    ├── run_boltz                         # Run Boltz-1
    └── run_chai                          # Run Chai-1
```


## Pipeline Steps

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

### 4. Energy Minimization (Not yet developed)

### 5. Extract Model Outputs

Extract predicted structures from model outputs:

```bash
bash ./scripts/extract_model_output.sh <dataset> <model_type>
```

### 6. Align Predicted Structures

Align predicted structures to the reference structures:

```bash
bash ./scripts/complex_structure_alignment.sh <dataset> <model_type>
```

### 7. Calculate Benchmark Result

Calculate evaluation metrics using PoseBusters:

```bash
bash ./scripts/calculate_benchmark_result.sh <dataset> <model_type>
```

## Environment Setup

### For Chai-1
```bash
pip install -r environments/chai.txt
```

### For Boltz-1
```bash
pip install -r environments/boltz.txt
```


