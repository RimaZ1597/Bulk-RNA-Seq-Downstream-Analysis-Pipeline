#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    COMBAT-SEQ BATCH CORRECTION WORKFLOW
========================================================================================
*/

// Import the ComBat-seq process
include { COMBATSEQ as COMBATSEQ_PROCESS } from '../modules/local/combat_seq/main.nf'

workflow COMBATSEQ {
    take:
    combatseq_input // tuple([meta], counts_file, metadata_file, config_file, experiment_id, outdir)
    
    main:
    // Execute ComBat-seq batch correction process
    COMBATSEQ_PROCESS(combatseq_input)
    
    emit:
    corrected_counts = COMBATSEQ_PROCESS.out.corrected_counts
    original_counts = COMBATSEQ_PROCESS.out.original_counts  
    metadata = COMBATSEQ_PROCESS.out.metadata
    log_files = COMBATSEQ_PROCESS.out.log_files
    combatseq_applied = COMBATSEQ_PROCESS.out.combatseq_applied  
}