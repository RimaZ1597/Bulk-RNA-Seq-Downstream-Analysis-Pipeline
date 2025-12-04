#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PATHWAY_ANALYSIS } from '../modules/local/pathway_analysis/main'

workflow {
    // Standalone pathway analysis workflow (CSV only, no GO)
    deseq2_results_dir = file(params.deseq2_results)
    config_yaml = file(params.config_yaml)
    
    ch_input = Channel.of([
        [id: 'pathway_analysis'],
        deseq2_results_dir,
        config_yaml
    ])
    
    PATHWAY_ANALYSIS(ch_input)
    
    // Output information
    PATHWAY_ANALYSIS.out.log.view { log_file ->
        "✅ Pathway analysis completed for study: ${params.study_name} (CSV only, no GO enrichment, all pathways)"
    }
}
