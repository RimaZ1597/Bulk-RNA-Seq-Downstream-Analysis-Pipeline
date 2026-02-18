#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    RNA-SEQ ANALYSIS PIPELINE
========================================================================================
*/

include { COMBATSEQ }      from './workflows/combatseq'
include { RNASEQPOWER }    from './workflows/rnaseqpower'
include { LIMMA_WORKFLOW } from './workflows/limma'
include { PATHWAY_WORKFLOW } from './workflows/pathway_analysis'

/*
========================================================================================
    PARAMETER DEFINITIONS
========================================================================================
*/

params.config_file  = "${projectDir}/interactive_setup/analysis_config.yaml"
params.project_name = null
params.outdir       = null
params.help         = false

def helpMessage() {
    log.info"""
    ================================================================================
                        RNA-seq Analysis Pipeline
    ================================================================================
    Usage:
      nextflow run main.nf --config_file <config.yaml> [options]

    Required Arguments:
      --config_file     Path to analysis configuration YAML file

    Optional Arguments:
      --project_name    Override project name from config
      --outdir          Override output directory from config
      --help            Show this help message

    Example:
      nextflow run main.nf \\
        --config_file interactive_setup/analysis_config.yaml \\
        --project_name PLX182687

    ================================================================================
    """.stripIndent()
}

def validateParams() {
    def errors = []
    if (!file(params.config_file).exists()) {
        errors << "Configuration file not found: ${params.config_file}"
    }
    if (errors) {
        log.error "Parameter validation failed:"
        errors.each { error -> log.error "  - ${error}" }
        log.error "\nRun with --help for usage information"
        System.exit(1)
    }
}
/* 
========================================================================================
    HELPER FUNCTIONS
========================================================================================
*/

def create_contrast_mapping(contrast_name) {
    """
    Create human-readable display name and safe file name from contrast.
    Returns: [original, display, safe_file_name]
    
    Example:
        Input:  "donor:19503_shear_stress:High (15 dynes/cm2)_vs_Static (0 dynes/cm2)"
        Output: 
            [0] original: "donor:19503_shear_stress:High (15 dynes/cm2)_vs_Static (0 dynes/cm2)"
            [1] display:  "donor 19503 shear_stress High (15 dynes/cm2) vs Static (0 dynes/cm2)"
            [2] safe:     "donor_19503_shear_stress_High_15_dynes_cm2_vs_Static_0_dynes_cm2"
    """
    def display = contrast_name.replaceAll(":", " ").replaceAll("_vs_", " vs ")
    
    def safe = contrast_name
        .replaceAll(/\([^)]*\)/, { it.replaceAll(/[()\/\s]/, '_') })
        .replaceAll(/[\/:,\s\(\)]/, '_')
        .replaceAll(/_+/, '_')
        .replaceAll(/^_|_$/, '')
    
    return [contrast_name, display, safe]
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    // Check help parameter first
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Validate parameters
    validateParams()
    
    // Parse configuration
    def yamlSlurper = new groovy.yaml.YamlSlurper()
    def cfg = yamlSlurper.parse(file(params.config_file))
    
    // Extract configuration
    def experiment_id = params.project_name ?: cfg.project_name ?: 'RNA_seq_Analysis'
    def outdir = params.outdir ?: cfg.analysis_parameters?.results_structure ?: 
                 "data/${experiment_id}/analysis_results"
    
    // Validate required configuration
    if (!cfg.input_files?.counts || !cfg.input_files?.metadata) {
        error """
        Configuration must define:
          - input_files.counts
          - input_files.metadata
        Please check your analysis_config.yaml file
        """.stripIndent()
    }
    
    // Read analysis parameters
    def perform_batch_correction    = cfg.analysis_parameters?.perform_batch_correction ?: false
    def perform_power_analysis      = cfg.analysis_parameters?.perform_power_analysis ?: false
    def perform_diff_expression     = cfg.analysis_parameters?.perform_differential_expression ?: false
    def perform_pathway_analysis    = cfg.analysis_parameters?.perform_pathway_analysis ?: false
    def de_tool                     = cfg.analysis_parameters?.de_tool ?: 'Not specified'
    
    // Read metadata column configuration
    def condition_col = cfg.experimental_design?.primary_factor ?: 
                        cfg.analysis_parameters?.metadata_columns?.condition ?: 'condition'
    def batch_col = cfg.experimental_design?.batch_factor ?: 
                    cfg.analysis_parameters?.metadata_columns?.batch
    
    // If batch correction enabled but no batch factor specified
    if (perform_batch_correction && !batch_col) {
        batch_col = 'group'
        log.info "ℹ️  Batch correction enabled: Using 'group' column for batch size validation"
    }
    
    // Read power analysis parameters
    def power_config = cfg.analysis_parameters?.power_analysis ?: [:]
    def effect_sizes = power_config.effect_sizes?.join(',') ?: '1.25,1.5,1.75,2'
    def power_levels = power_config.power_levels?.join(',') ?: '0.8,0.9'
    def alpha = power_config.alpha ?: 0.05
    def min_count = power_config.min_count ?: 10
    def min_prop = power_config.min_prop ?: 0.5
    
    // Pipeline header
    log.info """
    ================================================================================
                        RNA-SEQ ANALYSIS PIPELINE
    ================================================================================
    Experiment ID    : ${experiment_id}
    Output Directory : ${outdir}
    
    Input Files:
      Counts         : ${cfg.input_files.counts}
      Metadata       : ${cfg.input_files.metadata}
    
    Metadata Configuration:
      Condition column : ${condition_col}
      Batch column     : ${batch_col ?: 'None (global analysis)'}
    
    Differential Expression Tool:
      Selected tool    : ${de_tool}
      ${de_tool == 'DESeq2' ? '└─ Method: Negative binomial GLM with Wald test' : ''}
      ${de_tool == 'limma' ? '└─ Method: limma-voom with empirical Bayes moderation' : ''}
    
    Power Analysis Configuration:
      Effect sizes     : ${effect_sizes}
      Power levels     : ${power_levels}
      Alpha            : ${alpha}
      Min count        : ${min_count}
      Min proportion   : ${min_prop}
    
    Analysis Steps:
      ${perform_batch_correction ? '✓' : '⏸'} Step 1: Batch Correction (ComBat-seq)
      ${perform_power_analysis ? '✓' : '⏸'} Step 2: Power Analysis (RNASeqPower)
      ${perform_diff_expression ? '✓' : '⏸'} Step 3: Differential Expression (${de_tool})
      ${perform_pathway_analysis ? '✓' : '⏸'} Step 4: Pathway Analysis (FGSEA, MSigDB Collections & HumanGEM)
    ================================================================================
    """.stripIndent()
    
    // Define paths relative to projectDir
    def data_base_dir = "${projectDir}/interactive_setup/data"
    def project_data_dir = "${data_base_dir}/${experiment_id}"
    def counts_path = "${project_data_dir}/counts.csv"
    def metadata_path = "${project_data_dir}/metadata.csv"
    
    // Create channels with validation
    ch_counts = channel.fromPath(counts_path, checkIfExists: true)
    ch_metadata = channel.fromPath(metadata_path, checkIfExists: true)
    ch_config = channel.fromPath(params.config_file, checkIfExists: true)
    
    /*
     * ========================================================================
     * STEP 1: DATA PREPROCESSING (Filtering + Optional Batch Correction)
     * ========================================================================
     */
    log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    if (perform_batch_correction) {
        log.info "STARTING: Step 1 - Gene Filtering + Batch Correction (ComBat-seq)"
    } else {
        log.info "STARTING: Step 1 - Gene Filtering (Batch Correction Disabled)"
    }
    log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    // Always run the ComBat-seq module for filtering, regardless of batch correction setting
    // The module will handle batch correction based on the config file settings
    ch_combatseq_input = ch_counts
        .combine(ch_metadata)
        .combine(ch_config)
        .map { counts, metadata, config ->
            tuple(
                [
                    id: experiment_id,
                    step: 'data_preprocessing',
                    outdir: outdir
                ],
                counts,
                metadata,
                config,
                experiment_id,
                outdir
            )
        }
    
    COMBATSEQ(ch_combatseq_input)
    
    ch_corrected_counts  = COMBATSEQ.out.corrected_counts
    ch_original_counts   = COMBATSEQ.out.original_counts
    ch_metadata_out      = COMBATSEQ.out.metadata
    ch_combatseq_applied = COMBATSEQ.out.combatseq_applied
    
    if (perform_batch_correction) {
        COMBATSEQ.out.corrected_counts.first().view { "✓ Step 1 completed: Gene filtering + Batch correction" }
    } else {
        COMBATSEQ.out.corrected_counts.first().view { "✓ Step 1 completed: Gene filtering (batch correction skipped)" }
    }
    
    /*
     * ========================================================================
     * STEP 2: POWER ANALYSIS (RNASeqPower)
     * ========================================================================
     */
    if (perform_power_analysis) {
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log.info "STARTING: Step 2 - Power Analysis (RNASeqPower)"
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        if (perform_batch_correction) {
            log.info "  ℹ️  Using batch-corrected counts vs original filtered counts"
        } else {
            log.info "  ℹ️  Using filtered counts vs original filtered counts (no batch correction)"
        }
        
        ch_power_input = ch_original_counts
            .join(ch_corrected_counts)
            .join(ch_metadata_out)
            .map { meta, orig_file, corr_file, md_file ->
                tuple(
                    [
                        id: meta.id,
                        step: 'power_analysis',
                        outdir: outdir,
                        condition_column: condition_col,
                        batch_column: null,
                        effect_sizes: effect_sizes,
                        power_levels: power_levels,
                        alpha: alpha,
                        min_count: min_count,
                        min_prop: min_prop
                    ],
                    orig_file,
                    corr_file,
                    md_file
                )
            }
        
        RNASEQPOWER(ch_power_input)
        
        RNASEQPOWER.out.summary.first().view { meta, file ->
            "✓ Step 2 completed: Power analysis\n" +
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" +
            "📉 Power Analysis Summary Generated\n" +
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" +
            "Experiment : ${meta.id}\n" +
            "Summary    : ${file.name}\n" +
            "Location   : data/${meta.id}/analysis_results/power_analysis/\n" +
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        }
    } else {
        log.info "⏸  Step 2 skipped: Power analysis disabled"
    }
    
    /*
     * ========================================================================
     * STEP 3: DIFFERENTIAL EXPRESSION (limma-voom or DESeq2)
     * ========================================================================
     */
    if (perform_diff_expression) {
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        log.info "STARTING: Step 3 - Differential Expression (${de_tool})"
        log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        
        if (!cfg.comparisons || cfg.comparisons.size() == 0) {
            log.warn "⚠️  No contrasts defined - skipping differential expression"
            log.info "⏸  Step 3 skipped: No contrasts defined"
        } else {
            log.info "  Contrasts defined: ${cfg.comparisons.size()}"
            log.info "  DE tool: ${de_tool}"
            
            if (perform_batch_correction) {
                log.info "  ℹ️  Using batch-corrected counts for differential expression"
            } else {
                log.info "  ℹ️  Using filtered counts for differential expression (batch correction disabled)"
            }
            
            // Prepare input channel
            ch_de_input = ch_corrected_counts
                .join(ch_metadata_out)
                .join(ch_combatseq_applied)
                .combine(ch_config)
                .map { meta, counts_file, metadata_file, applied_flag_file, config_file ->
                    tuple(
                        meta.id,
                        counts_file,
                        metadata_file,
                        config_file,
                        applied_flag_file
                    )
                }
            
                        // Call appropriate workflow
            if (de_tool == 'limma') {
                LIMMA_WORKFLOW(ch_de_input)
                
                // ═══════════════════════════════════════════════════════
                // STEP 4: PATHWAY ANALYSIS (Triggered after LIMMA)
                // ═══════════════════════════════════════════════════════
                if (perform_pathway_analysis) {
                    log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                    log.info "STARTING: Step 4 - Pathway Analysis (FGSEA)"
                    log.info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                    
                    def subsystem_genes_path = "${projectDir}/subsystem_genes.csv"
                    
                    if (!file(subsystem_genes_path).exists()) {
                        log.error "❌ subsystem_genes.csv not found at: ${subsystem_genes_path}"
                        log.info "⏸  Step 4 skipped: Missing subsystem genes file"
                    } else {
                        // ✅ Create pathway input channel from limma results
                        ch_pathway_input = LIMMA_WORKFLOW.out.summary
                            .map { summary_file ->
                                // Add a small delay to ensure files are published
                                sleep(2000)
                                
                                // Use the published limma directory to find all contrast results  
                                def limma_base_dir = file("${projectDir}/data/${experiment_id}/analysis_results/limma")
                                
                                log.info "[PATHWAY] Scanning for contrast results in: ${limma_base_dir}"
                                
                                if (!limma_base_dir.exists()) {
                                    log.warn "[PATHWAY] Directory not found: ${limma_base_dir}"
                                    return []
                                }
                                
                                def pathway_inputs = []
                                limma_base_dir.eachDir { contrast_dir ->
                                    def all_genes_file = file("${contrast_dir}/all_genes_results.csv")
                                    if (all_genes_file.exists()) {
                                        def safe_contrast_name = contrast_dir.name
                                        def original_contrast = safe_contrast_name.replace("_", " ")
                                        
                                        log.info "[PATHWAY] Found contrast: ${safe_contrast_name} -> ${original_contrast}"
                                        
                                        pathway_inputs << [
                                            experiment_id,
                                            original_contrast,
                                            safe_contrast_name,
                                            all_genes_file,
                                            file(params.config_file),
                                            file(subsystem_genes_path)
                                        ]
                                    } else {
                                        log.warn "[PATHWAY] Missing file: ${all_genes_file}"
                                    }
                                }
                                
                                log.info "[PATHWAY] Created ${pathway_inputs.size()} pathway inputs"
                                log.info "[PATHWAY] Pathway inputs: ${pathway_inputs.collect { it[1] }.join(', ')}"
                                
                                return pathway_inputs
                            }
                            .flatMap { it }
                            .map { item -> 
                                log.info "[PATHWAY] ➤ Emitting individual contrast: ${item[1]} (${item[2]})"
                                tuple(item[0], item[1], item[2], item[3], item[4], item[5])
                            }
                        
                        PATHWAY_WORKFLOW(
                            ch_pathway_input
                        )
                        
                        PATHWAY_WORKFLOW.out.combined_summary.first().view { summary_file ->
                            "✓ Step 4 completed: Pathway Analysis (FGSEA)\n" +
                            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" +
                            "📊 Pathway Analysis Results Generated\n" +
                            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" +
                            "Experiment : ${experiment_id}\n" +
                            "Summary    : all_contrasts_pathway_summary.csv\n" +
                            "Location   : data/${experiment_id}/analysis_results/pathway_analysis/\n" +
                            "Contrasts  : ${cfg.comparisons.size()}\n" +
                            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                        }
                    }
                } else {
                    log.info "⏸  Step 4 skipped: Pathway Analysis disabled"
                }
                
                
                LIMMA_WORKFLOW.out.summary.first().view { summary_file ->
                    "✓ Step 3 completed: Differential Expression (limma-voom)\n" +
                    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" +
                    "📊 Differential Expression Results Generated\n" +
                    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n" +
                    "Experiment : ${experiment_id}\n" +
                    "Summary    : analysis_summary.csv\n" +
                    "Tool       : limma-voom\n" +
                    "Location   : data/${experiment_id}/analysis_results/limma/\n" +
                    "Contrasts  : ${cfg.comparisons.size()}\n" +
                    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                }
                
            } else if (de_tool == 'DESeq2') {
                log.warn "⚠️  DESeq2 workflow not yet implemented - coming soon!"
                log.info "⏸  Step 3 skipped: DESeq2 under development"
                
            } else {
                log.error "❌ Unknown DE tool: ${de_tool}"
                log.info "⏸  Step 3 skipped: Invalid DE tool"
            }
        }
    } else {
        log.info "⏸  Step 3 skipped: Differential Expression disabled"
    }


}

/*
========================================================================================
    WORKFLOW INTROSPECTION
========================================================================================
*/

workflow.onComplete {
    def yamlSlurper = new groovy.yaml.YamlSlurper()
    def cfg = yamlSlurper.parse(file(params.config_file ?: 'interactive_setup/analysis_config.yaml'))
    def exp_id = params.project_name ?: cfg.project_name ?: 'RNA_seq_Analysis'
    def de_tool = cfg.analysis_parameters?.de_tool ?: 'Not specified'
    def data_dir = "data/${exp_id}/analysis_results"
    
    println ""
    println "=" * 80
    println "PIPELINE COMPLETED"
    println "=" * 80
    println "Status           : ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}"
    println "Experiment ID    : ${exp_id}"
    println "DE Tool          : ${de_tool}"
    println "Duration         : ${workflow.duration}"
    println "Output Directory : ${data_dir}/"
    println ""
    
    if (workflow.success) {
        println "Results Summary:"
        if (file("${data_dir}/combatseq").exists()) {
            println "  • Batch Correction        : ${data_dir}/combatseq/"
        }
        if (file("${data_dir}/power_analysis").exists()) {
            println "  • Power Analysis          : ${data_dir}/power_analysis/"
        }
        if (file("${data_dir}/limma").exists()) {
            println "  • Differential Expression : ${data_dir}/limma/ (limma-voom)"
        } else if (file("${data_dir}/deseq2").exists()) {
            println "  • Differential Expression : ${data_dir}/deseq2/ (DESeq2)"
        }
        if (file("${data_dir}/pathway_analysis").exists()) {
            println "  • Pathway Analysis        : ${data_dir}/pathway_analysis/ (FGSEA)"
        }
        
        println ""
        println "Key Outputs:"
        if (file("${data_dir}/combatseq/combatseq_corrected_counts.csv").exists()) {
            println "  ✓ combatseq_corrected_counts.csv"
        }
        if (file("${data_dir}/power_analysis/power_comparison_summary.csv").exists()) {
            println "  ✓ power_comparison_summary.csv"
            println "  ✓ recommended_sample_sizes.csv"
        }
        if (file("${data_dir}/limma/analysis_summary.csv").exists()) {
            println "  ✓ analysis_summary.csv"
            println "  ✓ normalized_counts_cpm.csv"
            println "  ✓ all_contrasts_all_genes_combined.csv"
        }
        if (file("${data_dir}/deseq2/deseq2_summary.csv").exists()) {
            println "  ✓ deseq2_summary.csv (DESeq2 results)"
        }
        if (file("${data_dir}/pathway_analysis/all_contrasts_pathway_summary.csv").exists()) {
            println "  ✓ all_contrasts_pathway_summary.csv"
            println "  ✓ HumanGEM subsystem pathway results"
            println "  ✓ MSigDB pathway collections results"
        }
    } else {
        println "Error    : ${workflow.errorMessage ?: 'Unknown error'}"
        println "Work Dir : ${workflow.workDir}"
        println "Log file : ${workflow.workDir}/.nextflow.log"
    }
    println "=" * 80
    println ""
}

workflow.onError {
    println ""
    println "!" * 80
    println "WORKFLOW ERROR"
    println "!" * 80
    println "Message  : ${workflow.errorMessage}"
    println "Work Dir : ${workflow.workDir}"
    println "!" * 80
    println ""
}
