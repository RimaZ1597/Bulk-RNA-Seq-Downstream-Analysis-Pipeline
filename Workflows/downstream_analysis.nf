#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DESEQ2_ANALYSIS } from '../modules/local/deseq2/main'
include { PATHWAY_ANALYSIS } from '../modules/local/pathway_analysis/main'

workflow {
    // Combined downstream analysis workflow
    corrected_counts = file(params.corrected_counts)
    corrected_metadata = file(params.corrected_metadata) 
    config_yaml = file(params.config_yaml)
    
    // DESeq2 Analysis
    ch_deseq2_input = Channel.of([
        [id: 'deseq2_analysis'],
        corrected_counts,
        corrected_metadata,
        config_yaml
    ])
    
    DESEQ2_ANALYSIS(ch_deseq2_input)
    
    // Pathway Analysis
    ch_pathway_input = DESEQ2_ANALYSIS.out.deg_results
        .collect()
        .map { deseq2_files ->
            def deseq2_dir = deseq2_files[0].parent
            [
                [id: 'pathway_analysis'],
                deseq2_dir,
                file(params.config_yaml)
            ]
        }
    
    PATHWAY_ANALYSIS(ch_pathway_input)
    
    // Output information
    PATHWAY_ANALYSIS.out.log.view { log_file ->
        "Complete downstream analysis finished for study: ${params.study_name}"
    }
}
