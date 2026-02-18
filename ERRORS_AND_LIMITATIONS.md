# Pipeline Errors and Limitations - BulkRNAseq Limma Pipeline

## Overview
This document catalogs known errors, limitations, and troubleshooting information for the BulkRNAseq Pipeline based on execution logs and common issues encountered during development and testing.

## Current Known Errors

Based on analysis of actual execution logs (`.nextflow.log*` files), here are the documented errors:



### **ACTIVE ISSUE**: Incomplete Contrast Processing in Pathway Analysis
**Error Pattern**: Only one contrast folder generated in pathway analysis (Found in current executions)
```
[PATHWAY] Found contrast: treatment_GalXC_SLC25A5_vs_PBS_media_Fatty -> treatment GalXC SLC25A5 vs PBS media Fatty
[PATHWAY] Found contrast: treatment_PBS_vs_GalXC_SLC25A5_media_Lean -> treatment PBS vs GalXC SLC25A5 media Lean
[PATHWAY] ➤ Emitting individual contrast: treatment GalXC SLC25A5 vs PBS media Fatty (treatment_GalXC_SLC25A5_vs_PBS_media_Fatty)
# Only ONE pathway analysis task is submitted, missing second contrast
```

**Root Cause**: 
- Pipeline correctly identifies both contrasts (2 contrasts defined)
- Limma analysis generates folders for BOTH contrasts correctly
- Pathway analysis workflow only processes the FIRST contrast
- Second contrast (`treatment_PBS_vs_GalXC_SLC25A5_media_Lean`) is not submitted to pathway analysis

**Evidence from Current Execution**:
```
✓ Generated: limma/treatment_GalXC_SLC25A5_vs_PBS_media_Fatty/
✓ Generated: limma/treatment_PBS_vs_GalXC_SLC25A5_media_Lean/
✓ Generated: limma/treatment_GalXC_SLC25A5_vs_PBS_media_Fatty/pathway_analysis/
✗ Missing:   limma/treatment_PBS_vs_GalXC_SLC25A5_media_Lean/pathway_analysis/
```

**Debugging Steps**:
```bash
# Check if both contrasts are detected in limma output
ls -la data/PLX007704/analysis_results/limma/

# Should show BOTH:
# treatment_GalXC_SLC25A5_vs_PBS_media_Fatty/
# treatment_PBS_vs_GalXC_SLC25A5_media_Lean/

# Check pathway analysis workflow channel logic
grep -n "collect\|each\|flatten" workflows/pathway_analysis.nf

# Check if pathway analysis task was submitted for second contrast
grep "PATHWAY_ANALYSIS.*Lean" .nextflow.log
```

**Potential Root Causes**:
- Channel flattening issue in pathway workflow
- Contrast filtering logic excluding second contrast  
- Race condition in contrast directory scanning
- Channel collection not properly handling multiple contrasts

**Immediate Fix Needed**: 
- Investigate `workflows/pathway_analysis.nf` channel handling
- Ensure all contrasts from limma are passed to pathway analysis
- Verify channel collection and emission logic


### Memory Allocation Errors
**Error Pattern**: Out of memory errors during analysis
```
ERROR: Process LIMMA failed - Cannot allocate vector of size X GB
```
**Root Cause**:
- Insufficient memory allocation for large datasets
- Memory leaks in R processes

**Solution**:
```bash
# Increase memory limits
nextflow run main.nf --max_memory 128.GB

# Monitor memory usage
htop  # during execution
```

### 7. File Permission Issues
**Error Pattern**: Permission denied errors
```
ERROR: Permission denied - cannot write to output directory
```
**Root Cause**:
- Insufficient write permissions for output directories
- Network filesystem issues

**Solution**:
```bash
# Fix permissions
chmod -R 755 data/
chmod -R 755 results/

# Create directories manually if needed
mkdir -p data/PROJECT_ID/analysis_results
```

## Current Limitations

### 1. **ACTIVE LIMITATION**: Incomplete Multi-Contrast Pathway Analysis
**Current Issue**: Pipeline only processes one contrast for pathway analysis when multiple contrasts are defined

**Specific Problem**:
- Detects all contrasts correctly: "Contrasts defined: 2"  
- Limma analysis processes ALL contrasts successfully
- Pathway analysis only processes FIRST contrast
- Missing pathway results for subsequent contrasts

**Impact on Results**:
- Incomplete pathway enrichment analysis
- Missing biological insights for some contrasts
- Unbalanced results presentation
- Manual intervention required for complete analysis

**Status**: **NEEDS IMMEDIATE ATTENTION** - This affects current pipeline utility

### 2. Dataset Size Constraints
**Limitation**: Performance degradation with very large datasets
- **Genes**: Optimal performance up to ~30,000 genes
- **Samples**: Tested up to ~200 samples
- **Memory**: Requires 16-64 GB RAM depending on data size

**Impact**:
- Slower execution times for larger datasets
- Potential memory overflow errors
- Increased computational resource requirements

**Workarounds**:
```bash
# Use gene filtering to reduce dataset size
# Increase resource allocation
--max_memory 256.GB --max_cpus 32
```

### 2. Batch Effect Correction Limitations
**Limitation**: ComBat-seq requires specific experimental design constraints

**Requirements**:
- Minimum 2 samples per batch
- Conditions should not be completely confounded with batches
- Sufficient statistical power for batch effect estimation

**Common Issues**:
- Single-sample batches cause ComBat-seq failure
- Complete confounding leads to rank-deficient design matrix
- Insufficient replication prevents proper correction

**Current Behavior**:
- Pipeline gracefully falls back to filtered counts when ComBat-seq cannot be applied
- Warning messages indicate when batch correction is skipped

### 3. Pathway Analysis Database Dependencies
**Limitation**: Requires internet connectivity for pathway database access

**External Dependencies**:
- MSigDB collections (online access)
- Gene Ontology database updates
- KEGG pathway information

**Offline Limitations**:
- Cannot download latest pathway annotations
- May use cached/outdated pathway information
- HumanGEM analysis relies on local `subsystem_genes.csv`

### 4. Statistical Analysis Constraints
**Limitation**: Minimum sample size requirements for reliable statistics

**Requirements**:
- Minimum 3 samples per group for differential expression
- Power analysis requires sufficient dynamic range
- Pathway analysis needs adequate gene coverage

**Edge Cases**:
- Very small effect sizes may not be detectable
- Unbalanced designs reduce statistical power
- Missing pathway coverage for some gene sets



##  Planned Improvements

### 1. Error Handling Enhancements
- Better error messages with suggested solutions
- Automatic recovery for common failures
- Input validation before execution starts
- Resource requirement estimation

### 2. Performance Optimizations
- Streaming processing for large datasets
- Memory usage optimization
- Parallel processing improvements
- Caching strategies for repeated analyses

### 3. User Experience Improvements
- Interactive error diagnosis
- Progress indicators for long-running steps
- Automated troubleshooting suggestions
- Better logging and reporting
