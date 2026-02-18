# BulkRNAseq Pipeline - Bin Directory

This directory contains R scripts for bulk RNA-seq data analysis using the limma-voom framework. Each script is designed to work with the Nextflow pipeline and follows standardized input/output formats.

## Overview

The pipeline performs comprehensive RNA-seq analysis including:
- Batch effect correction with ComBat-seq
- Differential expression analysis with limma-voom
- Pathway enrichment analysis with FGSEA
- Power analysis for experimental design

## Scripts

### 1. combat_seq.R
**Purpose**: Batch effect correction using ComBat-seq

**Usage**:
```bash
Rscript combat_seq.R \
  --counts counts.csv \
  --metadata metadata.csv \
  --config analysis_config.yaml \
  --output_dir ./results
```

**Key Features**:
- Automatic batch effect detection and validation
- Gene filtering based on configuration settings
- Quality control metrics and visualization
- Graceful handling when batch correction is not needed
- Compatible with setup.py configuration format

**Inputs**:
- `counts.csv`: Raw count matrix (genes as rows, samples as columns)
- `metadata.csv`: Sample metadata with experimental factors
- `analysis_config.yaml`: Configuration from setup.py wizard

**Outputs**:
- `combatseq_corrected_counts.csv`: Batch-corrected counts
- `combatseq_original_counts.csv`: Filtered original counts
- `combatseq_metadata.csv`: Aligned metadata
- `combatseq_log.txt`: Detailed analysis log
- `combatseq_applied.txt`: Flag indicating if correction was applied

### 2. limma.R
**Purpose**: Differential expression analysis using limma-voom

**Usage**:
```bash
Rscript limma.R \
  --counts corrected_counts.csv \
  --metadata metadata.csv \
  --config analysis_config.yaml \
  --output_dir ./results \
  --contrast "specific_contrast_name"
```

**Key Features**:
- Multi-factor experimental design support
- Automatic contrast parsing from configuration
- TMM normalization with voom transformation
- Batch factor integration with blocking
- Comprehensive result files for each contrast
- Natural sorting for proper sample ordering

**Inputs**:
- `counts.csv`: Corrected count matrix from ComBat-seq
- `metadata.csv`: Sample metadata
- `analysis_config.yaml`: Configuration with contrast definitions

**Outputs**:
- `normalized_counts_cpm.csv`: TMM-normalized CPM values
- `voom_log2cpm.csv`: Voom-transformed expression values
- `analysis_summary.csv`: Summary statistics per contrast
- `all_genes_combined_results.csv`: Combined results across contrasts
- Contrast-specific folders with individual results

### 3. pathway_analysis.R
**Purpose**: Pathway enrichment analysis using FGSEA

**Usage**:
```bash
Rscript pathway_analysis.R \
  --limma_results all_genes_results.csv \
  --contrast_name "contrast_name" \
  --subsystem_genes subsystem_genes.csv \
  --config analysis_config.yaml \
  --output_dir ./results
```

**Key Features**:
- HumanGEM metabolic subsystem analysis
- MSigDB pathway collections (Hallmark, KEGG, GO, etc.)
- Ranked gene list creation using sign(logFC) × -log10(P-value)
- FDR correction and significance filtering
- Comprehensive pathway database coverage

**Inputs**:
- `all_genes_results.csv`: Limma results for specific contrast
- `subsystem_genes.csv`: HumanGEM subsystem gene mappings
- `analysis_config.yaml`: Project configuration

**Outputs**:
- `ranked_gene_list.csv`: Genes ranked by statistical significance
- `humangem_subsystem_pathways.csv`: HumanGEM pathway results
- `msigdb_*_pathways.csv`: MSigDB collection results
- `msigdb_all_collections_combined.csv`: Combined MSigDB results
- `pathway_analysis_summary.csv`: Summary statistics

### 4. power_analysis.R
**Purpose**: Statistical power analysis for experimental design

**Usage**:
```bash
Rscript power_analysis.R \
  --original_counts filtered_counts.csv \
  --corrected_counts corrected_counts.csv \
  --metadata metadata.csv \
  --condition_column condition \
  --batch_column batch \
  --outdir ./results
```

**Key Features**:
- RNASeqPower-based calculations
- Before/after batch correction comparison
- Sample size recommendations for future studies
- Effect size and power level analysis
- CV (coefficient of variation) assessment

**Inputs**:
- `original_counts.csv`: Pre-correction filtered counts
- `corrected_counts.csv`: Post-correction counts
- `metadata.csv`: Sample metadata with factors

**Outputs**:
- `power_distribution_original.csv`: Power analysis before correction
- `power_distribution_corrected.csv`: Power analysis after correction
- `power_comparison_original_vs_corrected.csv`: Improvement metrics
- `recommended_sample_sizes.csv`: Future study recommendations
- `power_analysis_log.txt`: Detailed interpretation guide

## Configuration Requirements

### analysis_config.yaml Structure
```yaml
project_name: "Your_Project"
experimental_design:
  primary_factor: "condition"
  batch_factor: "batch"  # Optional
  additional_factors: ["factor1", "factor2"]  # Optional
analysis_parameters:
  perform_differential_expression: true
  perform_batch_correction: true
comparisons:
  - "condition:treatment_vs_control"
  - "factor1:A_factor2:X_vs_Y"
filtering_config:
  method: "count_based"
  threshold: 10
  min_samples: 3
```

### Required Input Formats

#### Count Matrix (CSV)
```
gene,sample1,sample2,sample3
ENSG00000000003,1234,987,1456
ENSG00000000005,567,234,789
```

#### Metadata (CSV)
```
sample,condition,batch,other_factor
sample1,treatment,batch1,A
sample2,control,batch1,B
sample3,treatment,batch2,A
```

## Dependencies

### Required R Packages
```r
# Core analysis
library(edgeR)      # TMM normalization, DGEList
library(limma)      # limma-voom differential expression
library(sva)        # ComBat-seq batch correction

# Pathway analysis
library(fgsea)      # Fast gene set enrichment
library(msigdbr)    # MSigDB gene sets

# Power analysis
library(RNASeqPower) # Statistical power calculations

# Data handling
library(dplyr)      # Data manipulation
library(tidyr)      # Data reshaping
library(yaml)       # Configuration parsing
library(argparse)   # Command line arguments
```

## Pipeline Integration

These scripts are designed to work with the Nextflow pipeline:

1. **setup.py** → generates `analysis_config.yaml`
2. **combat_seq.R** → batch correction and filtering
3. **limma.R** → differential expression analysis
4. **pathway_analysis.R** → pathway enrichment (per contrast)
5. **power_analysis.R** → experimental design assessment

## Quality Control Features

### Automatic Validations
- Sample alignment between count matrix and metadata
- Natural sorting for proper sample ordering
- Missing value and infinite value handling
- Design matrix rank checking
- Batch effect feasibility assessment

### Error Handling
- Graceful degradation when batch correction fails
- Fallback to intercept-only models when needed
- Comprehensive logging for troubleshooting
- Input validation with informative error messages

## Output Organization

```
results/
├── combat_seq/
│   ├── combatseq_corrected_counts.csv
│   ├── combatseq_log.txt
│   └── combatseq_applied.txt
├── limma/
│   ├── normalized_counts_cpm.csv
│   ├── analysis_summary.csv
│   └── contrast_name/
│       ├── all_genes_results.csv
│       └── significant_genes.csv
├── pathway/
│   └── contrast_name/
│       ├── humangem_subsystem_pathways.csv
│       ├── msigdb_*.csv
│       └── pathway_analysis_summary.csv
└── power/
    ├── power_comparison_*.csv
    ├── recommended_sample_sizes.csv
    └── power_analysis_log.txt
```

## Best Practices

### Data Preparation
1. Ensure consistent gene identifiers (ENSEMBL/HGNC symbols)
2. Remove low-quality samples before analysis
3. Use appropriate filtering thresholds for gene expression
4. Validate experimental design balance

### Statistical Considerations
1. Minimum 3 samples per group for reliable statistics
2. Balance conditions across batches when possible
3. Consider effect sizes realistic for your biological system
4. Use appropriate multiple testing correction

### Troubleshooting
1. Check `*_log.txt` files for detailed diagnostics
2. Verify sample names match between files
3. Ensure sufficient samples per batch for ComBat-seq
4. Review design matrix for rank deficiency issues

## Version Information

- **Pipeline Version**: Compatible with setup.py configuration wizard
- **R Version**: Requires R ≥ 4.0
- **Key Dependencies**: limma ≥ 3.50, edgeR ≥ 3.36, sva ≥ 3.42
- **Last Updated**: December 2025

## Support

For issues or questions:
1. Check log files for specific error messages
2. Verify input file formats match specifications
3. Ensure all required R packages are installed
4. Review configuration file syntax

## License

This pipeline is part of the BulkRNAseq analysis framework developed for comprehensive RNA-seq data analysis in research environments.
