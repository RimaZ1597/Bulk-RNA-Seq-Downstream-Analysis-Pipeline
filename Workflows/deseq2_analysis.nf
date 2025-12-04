#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DESEQ2_ANALYSIS } from '../modules/local/deseq2/main'

workflow {
    // DESeq2 differential expression workflow (CSV only, contrast-based)
    corrected_counts = file(params.corrected_counts)
    metadata_file = file(params.metadata)
    config_file = file(params.config_yaml)
    
    ch_input = Channel.of([
        [id: 'deseq2_analysis'],
        corrected_counts,
        metadata_file,
        config_file
    ])
    
    DESEQ2_ANALYSIS(ch_input)
    
    // Output information
    DESEQ2_ANALYSIS.out.log.view { log_file ->
        "✅ DESeq2 analysis completed for study: ${params.study_name} (CSV outputs, organized by contrast)"
    }
}