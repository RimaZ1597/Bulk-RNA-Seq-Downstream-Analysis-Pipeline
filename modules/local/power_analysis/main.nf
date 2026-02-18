#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    POWER ANALYSIS MODULE
========================================================================================
    Performs RNA-seq power analysis using RNASeqPower Bioconductor package
    
    References:
    - Nextflow best practices: https://nextflow.io/docs/latest/process.html
    - RNASeqPower: https://bioconductor.org/packages/RNASeqPower/
========================================================================================
*/

process RNASEQPOWER {
    tag "${meta.id}"
    label 'process_medium'
    
    // Best practice: Organize outputs by type [1]
    publishDir "data/${meta.id}/analysis_results/power_analysis",
        mode: params.publish_dir_mode ?: 'copy',
        overwrite: true,
        pattern: "*.{csv,txt}",
        saveAs: { filename ->
            filename.endsWith('.txt') ? "logs/${filename}" : filename
        }
    
    // Best practice: Use container for reproducibility [1]
    container params.containers?.power_analysis ?: '/novo/projects/departments/compbio/sysbio/users/rvwx/containers/NF_BulkRNAseq.sif'
    
    input:
    tuple val(meta),
          path(original_counts, stageAs: 'original_counts.csv'),
          path(corrected_counts, stageAs: 'corrected_counts.csv'), 
          path(metadata, stageAs: 'metadata.csv')
    
    output:
    tuple val(meta), path("power_distribution_original.csv"),              emit: power_original
    tuple val(meta), path("power_distribution_corrected.csv"),             emit: power_corrected
    tuple val(meta), path("power_comparison_original_vs_corrected.csv"),   emit: power_comparison
    tuple val(meta), path("power_comparison_combined.csv"),                emit: power_combined
    tuple val(meta), path("recommended_sample_sizes.csv"),                 emit: recommendations
    tuple val(meta), path("power_comparison_summary.csv"),                 emit: summary
    tuple val(meta), path("power_analysis_log.txt"),                       emit: log_files
    path "versions.yml",                                                    emit: versions
    
    when:
    // Best practice: Allow conditional execution [2]
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    // Extract configuration - prioritize meta over task.ext (set in modules.config)
    def condition_col = meta.condition_column ?: task.ext.condition_column ?: 'condition'
    def batch_col = meta.batch_column ? "--batch_column ${meta.batch_column}" : (task.ext.batch_column ? "--batch_column ${task.ext.batch_column}" : "")
    def min_count = meta.min_count ?: task.ext.min_count ?: 10
    def min_prop = meta.min_prop ?: task.ext.min_prop ?: 0.5
    def effect_sizes = meta.effect_sizes ?: task.ext.effect_sizes ?: '1.25,1.5,1.75,2'
    def power_levels = meta.power_levels ?: task.ext.power_levels ?: '0.8,0.9'
    def alpha = meta.alpha ?: task.ext.alpha ?: 0.05
    
    """
    echo "═══════════════════════════════════════════════════════════════"
    echo "Starting Power Analysis: ${meta.id}"
    echo "═══════════════════════════════════════════════════════════════"
    echo "Configuration:"
    echo "  Condition column: ${condition_col}"
    echo "  Batch column: ${batch_col ?: 'None (global analysis)'}"
    echo "  Min count: ${min_count}"
    echo "  Min proportion: ${min_prop}"
    echo "  Effect sizes: ${effect_sizes}"
    echo "  Power levels: ${power_levels}"
    echo "  Alpha: ${alpha}"
    echo ""
    
    # Execute power analysis
    Rscript ${projectDir}/bin/power_analysis.R \\
        --original_counts ${original_counts} \\
        --corrected_counts ${corrected_counts} \\
        --metadata ${metadata} \\
        --experiment_id ${meta.id} \\
        --outdir . \\
        --condition_column ${condition_col} \\
        ${batch_col} \\
        --min_count ${min_count} \\
        --min_prop ${min_prop} \\
        --effect_sizes ${effect_sizes} \\
        --power_levels ${power_levels} \\
        --alpha ${alpha} \\
        ${args}
    
    # Best practice: Verify all expected outputs were created [2]
    echo ""
    echo "Verifying outputs..."
    
    for file in power_distribution_original.csv \\
                power_distribution_corrected.csv \\
                power_comparison_original_vs_corrected.csv \\
                power_comparison_combined.csv \\
                recommended_sample_sizes.csv \\
                power_comparison_summary.csv \\
                power_analysis_log.txt; do
        if [ ! -f "\$file" ]; then
            echo "ERROR: Expected output not created: \$file" >&2
            exit 1
        fi
        echo "  ✓ \$file"
    done
    
    # Best practice: Track software versions [1]
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-rnaseqpower: \$(Rscript -e "cat(as.character(packageVersion('RNASeqPower')))")
    END_VERSIONS
    
    echo ""
    echo "✅ Power analysis completed successfully"
    echo "═══════════════════════════════════════════════════════════════"
    """
    
    stub:
    // Best practice: Provide stub for testing [3]
    """
    # Create empty output files for testing
    touch power_distribution_original.csv
    touch power_distribution_corrected.csv
    touch power_comparison_original_vs_corrected.csv
    touch power_comparison_combined.csv
    touch recommended_sample_sizes.csv
    touch power_comparison_summary.csv
    echo "Stub execution for ${meta.id}" > power_analysis_log.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}