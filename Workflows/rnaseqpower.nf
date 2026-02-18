#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    RNASEQPOWER WORKFLOW
========================================================================================
*/

include { RNASEQPOWER as RNASEQPOWER_PROCESS } from '../modules/local/power_analysis/main.nf'

workflow RNASEQPOWER {
    
    take:
    ch_power_input
    
    main:
    
    // Validate inputs
    ch_power_input
        .map { meta, orig, corr, md ->
            if (!orig.exists()) {
                error "Original counts file not found: ${orig}"
            }
            if (!corr.exists()) {
                error "Corrected counts file not found: ${corr}"
            }
            if (!md.exists()) {
                error "Metadata file not found: ${md}"
            }
            if (orig.size() == 0) {
                error "Original counts file is empty: ${orig}"
            }
            if (corr.size() == 0) {
                error "Corrected counts file is empty: ${corr}"
            }
            if (md.size() == 0) {
                error "Metadata file is empty: ${md}"
            }
            tuple(meta, orig, corr, md)
        }
        .set { ch_validated }
    
    // Execute power analysis
    RNASEQPOWER_PROCESS(ch_validated)
    
    emit:
    power_original    = RNASEQPOWER_PROCESS.out.power_original
    power_corrected   = RNASEQPOWER_PROCESS.out.power_corrected
    power_comparison  = RNASEQPOWER_PROCESS.out.power_comparison
    power_combined    = RNASEQPOWER_PROCESS.out.power_combined
    recommendations   = RNASEQPOWER_PROCESS.out.recommendations
    summary           = RNASEQPOWER_PROCESS.out.summary
    log_files         = RNASEQPOWER_PROCESS.out.log
    versions          = RNASEQPOWER_PROCESS.out.versions
}