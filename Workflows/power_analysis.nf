#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { POWER_ANALYSIS } from '../modules/local/power_analysis/main'

workflow {
    // Standalone power analysis workflow (CSV only)
    original_counts = file(params.counts)
    corrected_counts = file(params.corrected_counts)
    corrected_metadata = file(params.corrected_metadata)
    config_yaml = file(params.config_yaml)
    
    ch_input = Channel.of([
        [id: 'power_analysis'],
        original_counts,
        corrected_counts,
        corrected_metadata,
        config_yaml
    ])
    
    POWER_ANALYSIS(ch_input)
    
    // Output information
    POWER_ANALYSIS.out.log.view { log_file ->
        "✅ Power analysis completed for study: ${params.study_name} (CSV outputs only)"
    }
}