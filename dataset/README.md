# PoseX


## Setup Environment

    conda env create -f environment.yml


## Generate Benchmark Dataset 
**Step 1:**  Download “Tabular Report - Entry IDs” that you are interested in from [RCSB PDB](https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_accession_info.initial_release_date%22%2C%22operator%22%3A%22greater_or_equal%22%2C%22negation%22%3Afalse%2C%22value%22%3A%222022-01-01%22%7D%7D%2C%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_accession_info.initial_release_date%22%2C%22operator%22%3A%22less_or_equal%22%2C%22negation%22%3Afalse%2C%22value%22%3A%222025-01-01%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_entry_info.selected_polymer_entity_types%22%2C%22operator%22%3A%22exact_match%22%2C%22negation%22%3Afalse%2C%22value%22%3A%22Protein%20(only)%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id%22%2C%22operator%22%3A%22exists%22%2C%22negation%22%3Afalse%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_entry_info.resolution_combined%22%2C%22operator%22%3A%22less_or_equal%22%2C%22negation%22%3Afalse%2C%22value%22%3A2%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%5D%2C%22label%22%3A%22text%22%7D%5D%7D%2C%22return_type%22%3A%22entry%22%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22results_content_type%22%3A%5B%22experimental%22%5D%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%2C%22scoring_strategy%22%3A%22combined%22%7D%2C%22request_info%22%3A%7B%22query_id%22%3A%2223a56d461e7e7e96f4065e59843158fe%22%7D%7D), you need to merge multiple *.txt files into one if the number of entries is greater than 10000. 
    
**Step 2:** Generate benchmark dataset 

    python main.py

*Inputs*:
- `--mode`: benchmark mode (self_dock or cross_dock).
- `--pdbid_path`: Path to the downloaded txt file containing Entry IDs in "Step 1".
- `--download_dir`: Folder to save the downloaded files.
- `--mmseqs_exec`: Path to the [MMseqs2](https://github.com/soedinglab/MMseqs2) binary (for protein clustering).

*Outputs*:
- dataset in a folder named `${mode}`


## Generate Custom Dataset 

    python create_custom_set.py

*Inputs*:
- `--name`: dataset name.
- `--pdb_ccd_path`: Path to the txt file containing items of f"{PDBID}_{CCDID}".
- `--download_dir`: Folder to save the downloaded files.

*Outputs*:
- dataset in a folder named `${name}`


## Directory Structure
```plaintext
.
├── create_custom_set.py    # Create custom dataset 
├── main.py                 # Entry point
├── posex
│   ├── align.py            # Cross alignment module
│   ├── ccd.py              # CCD utils module
│   ├── data.py             # Dataset generator module
│   ├── mmcif.py            # MMCIF parser module
│   ├── preprocess.py       # Dataset preprocessor module
│   └── utils.py            # Utils module
└── template
    ├── ccd_query.txt       # CCD query JSON template 
    ├── cross_dock.txt      # Cross-dock table Jinja2 template 
    ├── self_dock.txt       # Self-dock table Jinja2 template 
    └── vs_query.txt        # Validation score query JSON template 

```