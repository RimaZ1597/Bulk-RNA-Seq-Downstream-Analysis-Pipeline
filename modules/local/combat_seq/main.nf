#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    COMBAT-SEQ BATCH CORRECTION PROCESS
========================================================================================
    Following ComBat-seq GitHub repository implementation exactly
    
    Requirements:
    - Container: R image with sva, ComBat-seq packages
    - R script: bin/combat_seq.R (already exists)
    - Input: Raw counts + metadata
    - Output: Batch corrected counts and Filtered counts
========================================================================================
*/

process COMBATSEQ {
    tag "${meta.id}"
    label 'process_medium'
    
    publishDir "data/${meta.id}/analysis_results/combatseq", 
               mode: 'copy',
               overwrite: true,
               pattern: "*.{csv,txt}"

    input:
    tuple val(meta), path(counts), path(metadata), path(config_file), val(experiment_id), val(outdir)

    output:
    tuple val(meta), path("combatseq_corrected_counts.csv"), emit: corrected_counts
    tuple val(meta), path("combatseq_original_counts.csv"),  emit: original_counts
    tuple val(meta), path("combatseq_metadata.csv"),        emit: metadata
    tuple val(meta), path("combatseq_log.txt"),             emit: log_files
    tuple val(meta), path("combatseq_applied.txt"),         emit: combatseq_applied

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Execute ComBat-seq batch correction following GitHub implementation
    Rscript ${projectDir}/bin/combat_seq.R \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --config ${config_file} \\
        --output_dir . \\
        ${args}
    
    # Read the combatseq_applied flag written by the R script
    combatseq_flag=\$(cat combatseq_applied.txt)
    
    # Handle output files based on whether ComBat-seq was applied
    if [ "\$combatseq_flag" = "true" ]; then
        # ComBat-seq was applied - corrected counts file should exist
        if [ ! -f "combatseq_corrected_counts.csv" ]; then
            echo "Error: combatseq_corrected_counts.csv not generated when ComBat-seq was applied" >&2
            exit 1
        fi
    else
        # ComBat-seq was skipped - use filtered counts as the 'corrected' counts for pipeline consistency
        if [ ! -f "filtered_counts.csv" ]; then
            echo "Error: filtered_counts.csv not generated when ComBat-seq was skipped" >&2
            exit 1
        fi
        # Create symbolic link for pipeline compatibility
        ln -sf filtered_counts.csv combatseq_corrected_counts.csv
    fi
    
    if [ ! -f "combatseq_applied.txt" ]; then
        echo "Error: combatseq_applied.txt not generated" >&2
        exit 1
    fi
    """
    
    stub:
    """
    touch combatseq_corrected_counts.csv
    touch combatseq_original_counts.csv
    touch combatseq_metadata.csv
    touch combatseq_log.txt
    echo "false" > combatseq_applied.txt
    """
}