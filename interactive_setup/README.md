# Interactive Setup Directory - Configuration Wizard

## Purpose
**IMPORTANT**: This folder contains the interactive configuration wizard for setting up RNA-seq analysis parameters. The `data/` subfolder here is **SEPARATE** from the main pipeline `data/` folder and serves a different purpose.

## Folder Structure
```
interactive_setup/
├── setup.py              # Main configuration wizard
├── utils_io.py           # Input/output utilities  
├── utils_pluto.py        # Pluto.jl integration utilities
├── analysis_config.yaml  # Generated configuration file
├── config_summary.txt    # Human-readable config summary
├── pluto_token.txt       # Pluto.jl authentication token
├── data/                 # LOCAL data folder (for wizard testing)
│   ├── Experiment_ID/        # Sample datasets for configuration
│   ├── Experiment_ID/
│   └── *.csv             # Individual test files
└── results/              # Test analysis outputs
```

## Key Files

### setup.py - Configuration Wizard
**Purpose**: Interactive wizard that generates `analysis_config.yaml` for the pipeline
- **Usage**: `python setup.py`
- **Function**: Guides users through experimental design specification
- **Output**: Standardized configuration file for Nextflow pipeline

### data/ Subfolder vs Main data/
- **interactive_setup/data/**: Sample datasets for testing configuration wizard
- **../data/**: Main pipeline input data (separate, larger datasets)
- **Relationship**: The /data/folder serves as a testing ground for the configuration wizard to validate parameters before the main pipeline runs on the full datasets in the main ../data/ folder.



## Usage

```bash
# Run the configuration wizard
cd interactive_setup/
python setup.py

# Generated files are used by main pipeline
cd ..
nextflow run main.nf --config_file interactive_setup/analysis_config.yaml
```

## Output
- `analysis_config.yaml` → Used by main pipeline
- Configuration validates against sample data in local `data/` folder
- Main analysis runs on full datasets in `../data/` folder

## Status: **ACTIVE AND IMPORTANT**
This folder is essential for pipeline setup - it generates the configuration that drives the entire analysis workflow.
