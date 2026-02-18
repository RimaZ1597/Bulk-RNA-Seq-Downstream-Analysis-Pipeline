#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    LIMMA-VOOM DIFFERENTIAL EXPRESSION WORKFLOW  
========================================================================================
    Following Limma GitHub repository implementation exactly
    
    Input:  Batch-corrected counts + metadata + config + contrasts
    Output: Complete DEG tables (NO p-value filtering) + normalized counts
    
    Results Structure: {Experiment_ID}/Results/Differential_Expression/Limma/{Contrast}/
      ├── all_genes_results.csv          (Complete results - NO filtering)
      ├── voom_normalized_counts.csv     (Limma-voom normalized counts)  
      ├── limma_fit_object.rds          (R object for reproducibility)
      ├── limma_analysis_log.txt        (Analysis parameters and logs)
      └── contrast_summary.csv          (Summary statistics)
========================================================================================
*/

// workflow/limma.nf

nextflow.enable.dsl = 2

include { LIMMA } from '../modules/local/limma/main'

workflow LIMMA_WORKFLOW {
    take:
    combat_results_ch  // [project_id, counts.csv, metadata.csv, config.yaml, combatseq_applied]
    
    main:
    
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║         limma-voom Differential Expression Workflow       ║
    ╚═══════════════════════════════════════════════════════════╝
    """
    
    limma_input_ch = combat_results_ch.map { project_id, counts_csv, metadata_csv, config_yaml, combatseq_applied ->
        
        if (combatseq_applied == "true") {
            log.info "[${project_id}] Using ComBat-seq corrected counts"
        } else {
            log.info "[${project_id}] Using filtered counts (ComBat-seq skipped)"
        }
        
        tuple(project_id, counts_csv, metadata_csv, config_yaml, combatseq_applied)
    }
    
    // Add limma script as input
    limma_script_ch = channel.fromPath("${projectDir}/bin/limma.R", checkIfExists: true)
    
    LIMMA(limma_input_ch, limma_script_ch)
    
    emit:
    results = LIMMA.out.results
    summary = LIMMA.out.summary
    voom_plot = LIMMA.out.voom_plot
    normalized_log2cpm = LIMMA.out.normalized_log2cpm
    all_results = LIMMA.out.all_results
    log_file = LIMMA.out.log
}

workflow {
    if (!params.project_id) error "ERROR: --project_id required"
    if (!params.counts) error "ERROR: --counts required"
    if (!params.metadata) error "ERROR: --metadata required"
    if (!params.config) error "ERROR: --config required"
    
    input_ch = channel.of([
        params.project_id,
        file(params.counts),
        file(params.metadata),
        file(params.config),
        params.combatseq_applied ?: "unknown"
    ])
    
    LIMMA_WORKFLOW(input_ch)
}