# Workflows Directory - Nextflow DSL2 Workflows

## Overview
This directory contains the main workflow definitions for the BulkRNAseq pipeline using the limma-voom framework. Each workflow orchestrates specific analysis steps and coordinates the execution of modules for comprehensive RNA-seq data analysis.

## Workflow Structure
```
workflows/
├── combatseq.nf         # Batch effect correction workflow
├── limma.nf            # Differential expression analysis workflow  
├── pathway_analysis.nf  # Pathway enrichment analysis workflow
└── rnaseqpower.nf      # Statistical power analysis workflow
```

## Active Workflows

### 1. COMBATSEQ Workflow (`combatseq.nf`)
**Purpose**: Orchestrates batch effect correction using ComBat-seq

**Process Flow**:
1. Input validation (counts, metadata, configuration)
2. ComBat-seq batch correction execution
3. Quality control assessment
4. Output organization and publishing

**Key Features**:
- Automatic batch effect detection
- Gene filtering based on configuration
- Graceful handling when correction not needed
- Quality metrics and logging

**Input Channels**:
```groovy
ch_combatseq_input = [
    meta: [id: project_id, step: 'combat_seq'],
    counts: path(counts.csv),
    metadata: path(metadata.csv),
    config: path(analysis_config.yaml)
]
```

**Output Channels**:
```groovy
corrected_counts    // Batch-corrected count matrix
original_counts     // Filtered original counts  
metadata           // Aligned sample metadata
log_files          // Analysis logs and reports
applied_flag       // Boolean flag indicating if correction was applied
```

### 2. LIMMA Workflow (`limma.nf`)
**Purpose**: Coordinates differential expression analysis using limma-voom

**Process Flow**:
1. Input preparation and validation
2. TMM normalization and voom transformation
3. Design matrix creation and contrast fitting
4. Statistical analysis and result generation
5. Output organization by contrast

**Key Features**:
- Multi-factor experimental design support
- Automatic contrast parsing from configuration
- Batch factor integration with blocking
- Comprehensive result files per contrast
- Natural sample ordering

**Input Channels**:
```groovy
ch_limma_input = [
    meta: [id: project_id, step: 'limma'],
    corrected_counts: path(corrected_counts.csv),
    metadata: path(metadata.csv),
    config: path(analysis_config.yaml),
    combatseq_flag: path(combatseq_applied.txt)
]
```

**Output Channels**:
```groovy
normalized_counts   // TMM-normalized and voom-transformed data
analysis_summary   // Summary statistics per contrast
combined_results   // All contrasts combined results
contrast_folders   // Individual contrast result directories
r_objects          // Saved R objects for reproducibility
```

### 3. PATHWAY_ANALYSIS Workflow (`pathway_analysis.nf`)
**Purpose**: Manages pathway enrichment analysis using FGSEA

**Process Flow**:
1. Limma results processing and gene ranking
2. HumanGEM subsystem pathway analysis
3. MSigDB pathway collection analysis
4. Result compilation and summary generation

**Key Features**:
- Multiple pathway database support
- Per-contrast analysis execution
- FGSEA-based enrichment with FDR correction
- Comprehensive pathway coverage

**Input Channels**:
```groovy
ch_pathway_input = [
    meta: [id: project_id, contrast: contrast_name],
    limma_results: path(all_genes_results.csv),
    subsystem_genes: path(subsystem_genes.csv),
    config: path(analysis_config.yaml)
]
```

**Output Channels**:
```groovy
ranked_genes       // Genes ranked by statistical significance
pathway_results    // Enrichment results per database
summary_stats      // Analysis summary and statistics
```

### 4. RNASEQPOWER Workflow (`rnaseqpower.nf`)
**Purpose**: Conducts statistical power analysis for experimental design

**Process Flow**:
1. Input validation and preparation
2. Power analysis on original counts
3. Power analysis on corrected counts
4. Comparison and improvement assessment
5. Sample size recommendations

**Key Features**:
- Before/after batch correction comparison
- Multiple effect size and power level testing
- Sample size recommendations for future studies
- CV (coefficient of variation) assessment

**Input Channels**:
```groovy
ch_power_input = [
    meta: [id: project_id, step: 'power_analysis'],
    original_counts: path(original_counts.csv),
    corrected_counts: path(corrected_counts.csv),
    metadata: path(metadata.csv),
    config: path(analysis_config.yaml)
]
```

**Output Channels**:
```groovy
power_original     // Power analysis before correction
power_corrected    // Power analysis after correction
power_comparison   // Improvement metrics
power_combined     // Combined analysis results
recommendations    // Sample size recommendations
summary           // Analysis summary statistics
log_files         // Detailed logs and interpretation
```

## Workflow Integration Patterns

### Sequential Execution
The main pipeline orchestrates workflows in sequence:
```groovy
// Main pipeline flow
COMBATSEQ(ch_input)
LIMMA(COMBATSEQ.out.corrected)
PATHWAY_ANALYSIS(LIMMA.out.contrasts)
RNASEQPOWER(COMBATSEQ.out.original, COMBATSEQ.out.corrected)
```

### Parallel Processing
Some workflows support parallel execution:
```groovy
// Pathway analysis per contrast (parallel)
LIMMA.out.contrasts
    .flatten()
    .map { contrast -> [meta + [contrast: contrast.name], contrast] }
    .set { ch_pathway_per_contrast }

PATHWAY_ANALYSIS(ch_pathway_per_contrast)
```

### Conditional Execution
Workflows can be conditionally executed based on configuration:
```groovy
// Conditional workflow execution
if (perform_batch_correction) {
    COMBATSEQ(ch_input)
} else {
    // Skip to limma with original counts
}

if (perform_power_analysis) {
    RNASEQPOWER(ch_power_input)
}
```

## Channel Management

### Standard Channel Structure
All workflows follow consistent channel patterns:
```groovy
// Input channel pattern
tuple val(meta), path(input_files), path(config)

// Output channel pattern
tuple val(meta), path(output_files), emit: channel_name
```

### Metadata Propagation
Metadata is consistently propagated through workflows:
```groovy
meta = [
    id: "Experiment_ID",           // Project identifier
    step: "limma",             // Current analysis step
    contrast: "High_vs_Static", // Specific contrast (if applicable)
    timestamp: timestamp        // Processing timestamp
]
```

### Error Handling
Workflows include robust error handling:
```groovy
// Input validation
ch_input
    .map { meta, files ->
        files.each { file ->
            if (!file.exists()) {
                error "Required file not found: ${file}"
            }
        }
        return [meta, files]
    }
```

## Configuration Integration

### Analysis Configuration
Workflows read parameters from `analysis_config.yaml`:
```yaml
analysis_parameters:
  perform_batch_correction: true
  perform_differential_expression: true
  perform_pathway_analysis: true
  perform_power_analysis: true

experimental_design:
  primary_factor: "condition"
  batch_factor: "batch"
  additional_factors: ["donor", "plate"]

comparisons:
  - "condition:treatment_vs_control"
  - "donor:A_condition:high_vs_low"
```

### Workflow-Specific Parameters
Each workflow can access specific configuration sections:
```groovy
// In limma.nf
def limma_config = config.analysis_parameters.limma ?: [:]
def fdr_threshold = limma_config.fdr_threshold ?: 0.05

// In pathway_analysis.nf  
def pathway_config = config.analysis_parameters.pathway_analysis ?: [:]
def min_pathway_size = pathway_config.min_size ?: 15
```

## Resource Management

### Workflow-Level Resources
Workflows can specify resource requirements:
```groovy
// In workflow definition
workflow LIMMA {
    // Resource hints for workflow scheduling
    maxForks = 4
    errorStrategy = 'retry'
    maxRetries = 2
}
```

### Module Resource Allocation
Resources are primarily managed at the module level via `conf/modules.config`:
```groovy
process {
    withName: 'COMBATSEQ' {
        cpus = 2
        memory = '12.GB'
        time = '6.h'
    }
    
    withName: 'LIMMA' {
        cpus = 4
        memory = '16.GB'
        time = '8.h'
    }
}
```

## Output Organization

### Standardized Output Structure
All workflows follow consistent output organization:
```
data/{project_id}/analysis_results/
├── combatseq/              # Batch correction results
│   ├── corrected_counts.csv
│   └── log.txt
├── limma/                  # Differential expression results
│   ├── analysis_summary.csv
│   ├── normalized_counts.csv
│   └── {contrast_name}/
│       ├── all_genes_results.csv
│       └── significant_genes.csv
├── pathway_analysis/       # Pathway enrichment results
│   └── {contrast_name}/
│       ├── humangem_pathways.csv
│       └── msigdb_pathways.csv
└── power_analysis/         # Power analysis results
    ├── power_comparison.csv
    └── recommendations.csv
```

### Publishing Strategy
Workflows use consistent publishing patterns:
```groovy
publishDir "data/${meta.id}/analysis_results/${meta.step}", 
           mode: 'copy',
           pattern: "*.{csv,txt,png,pdf}"
```

## Testing and Validation

### Workflow Testing
Each workflow can be tested independently:
```bash
# Test individual workflow
nextflow run workflows/combatseq.nf \
  --input test_data/test_input.csv \
  --config test_data/test_config.yaml

# Test with minimal data
nextflow run workflows/limma.nf -profile test
```

### Integration Testing
Full workflow integration testing:
```bash
# Test complete pipeline
nextflow run main.nf -profile test,docker \
  --config_file test_configs/minimal_config.yaml
```

### Validation Checks
Workflows include built-in validation:
- Input file existence and format checking
- Configuration parameter validation
- Output file generation verification
- Statistical result sanity checks

## Development Guidelines

### Adding New Workflows
1. **Create Workflow File**: Follow naming convention (`analysis_type.nf`)
2. **Define Interface**: Standardize input/output channels
3. **Include Modules**: Import required modules from `modules/local/`
4. **Add Configuration**: Update `conf/` files as needed
5. **Document**: Update this README and add inline documentation
6. **Test**: Create test cases and validate functionality

### Workflow Template
```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import required modules
include { MODULE_NAME } from '../modules/local/module_name/main.nf'

workflow WORKFLOW_NAME {
    take:
    ch_input
    
    main:
    // Input validation
    ch_validated = ch_input.map { meta, files ->
        // Validation logic
        return [meta, files]
    }
    
    // Execute modules
    MODULE_NAME(ch_validated)
    
    emit:
    results = MODULE_NAME.out.results
    log     = MODULE_NAME.out.log
}
```

### Best Practices
1. **Consistent Naming**: Use descriptive workflow and channel names
2. **Error Handling**: Include comprehensive input validation
3. **Documentation**: Add clear comments and docstrings
4. **Modularity**: Keep workflows focused and modular
5. **Testing**: Include test cases for edge conditions

## Troubleshooting

### Common Issues

#### 1. Channel Formatting Errors
```bash
# Check channel structure
.view { it -> "Channel content: ${it}" }
```

#### 2. Module Import Problems
```bash
# Verify module paths
include { MODULE } from '../modules/local/module/main.nf'
```

#### 3. Configuration Access Issues
```bash
# Debug configuration loading
config.view()
```

#### 4. Resource Allocation Problems
```bash
# Check resource requirements in modules.config
# Increase memory/CPU allocation if needed
```

### Debugging Workflows
```bash
# Run with debug output
nextflow run workflow.nf -with-trace -with-report

# Check specific workflow execution
nextflow log -f name,status,duration,realtime

# Resume failed execution
nextflow run workflow.nf -resume
```

## Performance Optimization

### Workflow Optimization
- **Parallel Execution**: Identify parallelizable steps
- **Resource Tuning**: Optimize CPU/memory allocation
- **Caching Strategy**: Use Nextflow work directory caching
- **Container Optimization**: Use efficient container images

### Monitoring
- **Execution Reports**: Generate timeline and trace reports
- **Resource Usage**: Monitor CPU and memory consumption
- **Bottleneck Identification**: Identify slow steps for optimization

## Version Information

- **Nextflow DSL**: DSL2
- **Nextflow Version**: ≥22.10.1
- **Workflow Standard**: nf-core compatible
- **Last Updated**: December 2025

## Support

For workflow-related issues:
1. **Execution Problems**: Check `.nextflow.log` and process logs
2. **Channel Issues**: Use `.view()` to debug channel content
3. **Module Integration**: Verify module interfaces and imports
4. **Resource Issues**: Adjust allocations in `conf/modules.config`
5. **Configuration Problems**: Validate YAML syntax and parameter values

Each workflow is designed to be modular, testable, and maintainable while following nf-core best practices for reproducible bioinformatics analysis.
