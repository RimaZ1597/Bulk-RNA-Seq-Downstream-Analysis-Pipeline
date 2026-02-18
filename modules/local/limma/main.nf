// modules/local/limma/main.nf

nextflow.enable.dsl = 2

/*
========================================================================================
    LIMMA-VOOM DIFFERENTIAL EXPRESSION MODULE
========================================================================================
    Following Limma GitHub implementation exactly
    
    Requirements:
    - Container: Bioconductor image with limma, edgeR packages
    - R script: bin/limma.R (standardized interface)
    - Input: Batch-corrected counts + metadata + contrasts
    - Output: Complete DEG tables (NO p-value filtering) + normalized counts
========================================================================================
*/

process LIMMA {
    tag "${project_id}"
    label 'process_medium'
    
    publishDir "data/${project_id}/analysis_results/limma", mode: 'copy'
    
    input:
    tuple val(project_id), path(counts_csv), path(metadata_csv), path(config_yaml), path(combatseq_applied_file)
    path limma_script
    
    output:
    path "**", emit: results
    path "analysis_summary.csv", emit: summary
    path "voom_mean_variance_plot.png", emit: voom_plot, optional: true
    path "normalized_counts_log2cpm.csv", emit: normalized_log2cpm
    path "all_contrasts_combined_results.csv", emit: all_results
    path "limma_log.txt", emit: log
    
    script:
    """
    # Read combatseq_applied flag from file
    combatseq_applied=\$(cat ${combatseq_applied_file})
    
    # Run limma analysis (publishDir handles the path structure)
    Rscript ${limma_script} \\
        --counts ${counts_csv} \\
        --metadata ${metadata_csv} \\
        --config ${config_yaml} \\
        --output_dir . \\
        --contrast all \\
        --combatseq_applied \$combatseq_applied \\
        2>&1 | tee limma_log.txt
    
    # Validate completion
    if [ \$? -ne 0 ]; then
        echo "ERROR: limma analysis failed" >&2
        exit 1
    fi
    
    if [ ! -f analysis_summary.csv ]; then
        echo "ERROR: analysis_summary.csv not created" >&2
        exit 1
    fi
    
    echo "✓ limma analysis completed successfully"
    """
}