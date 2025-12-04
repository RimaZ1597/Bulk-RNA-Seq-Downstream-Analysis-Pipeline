#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { COMBAT_SEQ } from '../modules/local/combat_seq/main'

workflow {
    // Standalone ComBat-seq workflow with study support
    counts_file = file(params.counts)
    metadata_file = file(params.metadata)
    config_file = file(params.config_yaml)
    
    ch_input = Channel.of([
        [id: 'batch_correction'], 
        counts_file, 
        metadata_file, 
        config_file
    ])
    
    COMBAT_SEQ(ch_input)
    
    // Output information
    COMBAT_SEQ.out.log.view { log_file ->
        "ComBat-seq batch correction completed for study: ${params.study_name}"
    }
}
