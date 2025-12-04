#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import all modules (adjust based on your actual module structure)
include { COMBAT_SEQ } from './modules/local/combat_seq/main'
include { POWER_ANALYSIS } from './modules/local/power_analysis/main'
include { DESEQ2_ANALYSIS } from './modules/local/deseq2/main'
include { PATHWAY_ANALYSIS } from './modules/local/pathway_analysis/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParam(param_name, param_value) {
    if (!param_value) {
        error "Parameter ${param_name} is required"
    }
    if (!file(param_value).exists()) {
        error "File specified by ${param_name} does not exist: ${param_value}"
    }
    return file(param_value)
}

def setDatasetParameters() {
    // Handle dataset-specific parameters
    def resolvedParams = [:]
    
    // Set defaults
    resolvedParams.counts = params.counts
    resolvedParams.metadata = params.metadata
    resolvedParams.config_yaml = params.config_yaml ?: 'conf/analysis_config.yaml'
    resolvedParams.study_name = params.study_name ?: 'default_study'
    
    // Check for dataset flags and override if not explicitly set
    if (params.LAMTOR_KD) {
        if (!params.counts) resolvedParams.counts = 'data/LAMTOR_KD/raw_counts.csv'
        if (!params.metadata) resolvedParams.metadata = 'data/LAMTOR_KD/metadata.csv'
        if (!params.config_yaml || params.config_yaml == 'conf/analysis_config.yaml') {
            resolvedParams.config_yaml = 'conf/analysis_config_LAMTOR_KD.yaml'
        }
        if (!params.study_name || params.study_name == 'default_study') {
            resolvedParams.study_name = 'LAMTOR_KD'
        }
        log.info "🎯 Using LAMTOR_KD dataset"
    } else if (params.REPSOX) {
        if (!params.counts) resolvedParams.counts = 'data/REPSOX/raw_counts.csv'
        if (!params.metadata) resolvedParams.metadata = 'data/REPSOX/metadata.csv'
        if (!params.config_yaml || params.config_yaml == 'conf/analysis_config.yaml') {
            resolvedParams.config_yaml = 'conf/analysis_config_REPSOX.yaml'
        }
        if (!params.study_name || params.study_name == 'default_study') {
            resolvedParams.study_name = 'REPSOX'
        }
        log.info "🎯 Using REPSOX dataset"
    } else if (params.Adipose_Endothelial_Cell_Flow) {
        if (!params.counts) resolvedParams.counts = 'data/Adipose_Endothelial_Cell_Flow/raw_counts.csv'
        if (!params.metadata) resolvedParams.metadata = 'data/Adipose_Endothelial_Cell_Flow/metadata.csv'
        if (!params.config_yaml || params.config_yaml == 'conf/analysis_config.yaml') {
            resolvedParams.config_yaml = 'conf/analysis_config_Adipose_Endothelial_Cell_Flow.yaml'
        }
        if (!params.study_name || params.study_name == 'default_study') {
            resolvedParams.study_name = 'Adipose_Endothelial_Cell_Flow'
        }
        log.info "🎯 Using Adipose_Endothelial_Cell_Flow dataset"
    }
    
    return resolvedParams
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW - COMPLETE RNA-SEQ ANALYSIS PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Handle dataset-specific parameters first
    def resolvedParams = setDatasetParameters()
    
    // Update global params with resolved values
    params.study_name = resolvedParams.study_name
    params.counts = resolvedParams.counts
    params.metadata = resolvedParams.metadata
    params.config_yaml = resolvedParams.config_yaml
    
    // Print parameter summary
    log.info """
    ===============================================
    🧬 RNA-SEQ DOWNSTREAM ANALYSIS PIPELINE
    ===============================================
    Counts matrix    : ${resolvedParams.counts ?: 'Not specified'}
    Metadata         : ${resolvedParams.metadata ?: 'Not specified'}  
    Configuration    : ${resolvedParams.config_yaml ?: 'conf/analysis_config.yaml'}
    Output directory : ${params.output_dir ?: 'results'}
    Study name       : ${resolvedParams.study_name ?: 'default_study'}
    
    Pipeline Steps:
    1. ComBat-seq Batch Correction
    2. Power Analysis (Enhanced)
    3. DESeq2 Differential Expression
    4. Pathway Analysis (Multi-Collection)
    ===============================================
    """
    
    // Input validation
    counts_file = checkPathParam('counts', resolvedParams.counts)
    metadata_file = checkPathParam('metadata', resolvedParams.metadata)
    config_file = checkPathParam('config_yaml', resolvedParams.config_yaml)
    
    // Create study-specific output directory
    study_output_dir = "${params.output_dir}/${resolvedParams.study_name}"
    log.info "Using counts file: ${counts_file}"
    log.info "Using metadata file: ${metadata_file}"
    log.info "Using config file: ${config_file}"
    log.info "Results will be saved to: ${study_output_dir}"

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 1: COMBAT-SEQ BATCH CORRECTION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    log.info "🔄 Starting Step 1: ComBat-seq Batch Correction..."
    
    ch_combat_input = Channel.of([
        [id: 'combat_correction'],
        counts_file,
        metadata_file,
        config_file
    ])
    
    COMBAT_SEQ(ch_combat_input)
    
    COMBAT_SEQ.out.corrected_data.view { meta, corrected_counts, corrected_metadata ->
        "✅ Step 1 Complete - ComBat-seq batch correction: ${corrected_counts}"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 2: POWER ANALYSIS (ENHANCED CSV OUTPUT)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    log.info "🔄 Starting Step 2: Enhanced Power Analysis (CSV outputs only)..."
    
    ch_power_input = Channel.of([
        [id: 'power_analysis'],
        counts_file,  // Original counts
    ]).join(
        COMBAT_SEQ.out.corrected_data.map { meta, corrected_counts, corrected_metadata ->
            [meta, corrected_counts, corrected_metadata]
        }, by: 0
    ).combine(Channel.of(config_file))
    .map { meta, original_counts, corrected_counts, corrected_metadata, config ->
        [
            meta,
            original_counts,
            corrected_counts,
            corrected_metadata,
            config
        ]
    }
    
    POWER_ANALYSIS(ch_power_input)
    
    POWER_ANALYSIS.out.log.view { log_file ->
        "✅ Step 2 Complete - Enhanced Power Analysis: ${log_file}"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 3: DIFFERENTIAL EXPRESSION ANALYSIS (BY CONTRAST)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    log.info "🔄 Starting Step 3: DESeq2 Differential Expression Analysis..."
    
    ch_deseq2_input = COMBAT_SEQ.out.corrected_data
        .combine(Channel.of(config_file))
        .map { meta, corrected_counts, corrected_metadata, config ->
            [
                meta,
                corrected_counts,
                corrected_metadata,
                config
            ]
        }
    
    DESEQ2_ANALYSIS(ch_deseq2_input)
    
    DESEQ2_ANALYSIS.out.log.view { log_file ->
        "✅ Step 3 Complete - DESeq2 Analysis: ${log_file}"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STEP 4: PATHWAY ANALYSIS (MULTI-COLLECTION)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    log.info "🔄 Starting Step 4: Enhanced Pathway Analysis..."
    
    // Fixed channel reference - use deg_results instead of deseq2_results
    ch_pathway_input = DESEQ2_ANALYSIS.out.deseq2_results
    
    PATHWAY_ANALYSIS(ch_pathway_input)
    
    PATHWAY_ANALYSIS.out.log.view { log_file ->
        "✅ Step 4 Complete - Enhanced Pathway Analysis: ${log_file}"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FINAL SUMMARY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_final_summary = COMBAT_SEQ.out.log
        .combine(POWER_ANALYSIS.out.log)
        .combine(DESEQ2_ANALYSIS.out.log)
        .combine(PATHWAY_ANALYSIS.out.log)
        .view { combat_log, power_log, deseq2_log, pathway_log ->
            """
            🎉 MULTI-DATASET RNA-SEQ PIPELINE COMPLETED SUCCESSFULLY! 🎉
            =============================================================
            📊 Results Summary for Study: ${params.study_name}
            
            ├── 01_Batch_correction/                    ✅ ComBat-seq correction completed
            ├── 02_Power_analysis/                      ✅ Enhanced power analysis completed
            │   ├── Combined dataset analysis
            │   ├── Sample size scenarios
            │   ├── Gene-level power calculations
            │   └── Individual vs combined comparisons
            └── 03_Differential_Expression_&_Pathway_Analysis/
                ├── [Contrast_Name]/
                │   ├── Differential_Expression_Analysis/  ✅ DESeq2 results (CSV only)
                │   └── Pathway_Analysis/                  ✅ Multi-collection pathway results
                │       ├── Hallmark pathways
                │       ├── KEGG pathways
                │       ├── Reactome pathways
                │       ├── GO pathways
                │       └── Regulatory targets
            
            📁 Check your results in: ${study_output_dir}/
            📄 Final reports: pathway_analysis_final_report.txt
            🎯 Output Format: CSV files only (no visualizations)
            🔬 Enhanced Features: Combined dataset power analysis, multi-collection pathways
            """
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW COMPLETION HOOKS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    // Get the resolved study name
    def resolvedParams = setDatasetParameters()
    
    if (workflow.success) {
        log.info """
        🎯 Multi-Dataset Pipeline execution completed successfully!
        📋 Execution Summary:
        - Study: ${resolvedParams.study_name}
        - Duration: ${workflow.duration}
        - Success: ${workflow.success}
        - Results directory: ${params.output_dir}/${resolvedParams.study_name}
        
        📊 Analysis Steps Completed:
        ✅ Step 1: ComBat-seq batch correction
        ✅ Step 2: Enhanced power analysis (combined dataset + scenarios)
        ✅ Step 3: DESeq2 differential expression (contrast-based)
        ✅ Step 4: Multi-collection pathway analysis
        
        🔬 Enhanced Features Delivered:
        - Combined dataset power analysis
        - Sample size scenario modeling
        - Gene-level power calculations
        - Multiple pathway collections (Hallmark, KEGG, Reactome, GO, etc.)
        - Contrast-specific organization
        - Comprehensive CSV-only outputs
        
        🎉 Complete enhanced analysis pipeline finished successfully!
        """.stripIndent()
    } else {
        log.error """
        ❌ Pipeline execution failed!
        Error details:
        - Exit status: ${workflow.exitStatus}
        - Error message: ${workflow.errorMessage}
        - Study: ${resolvedParams.study_name}
        """.stripIndent()
    }
}

workflow.onError {
    log.error "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}