process POWER_ANALYSIS {
    tag "Power analysis (CSV only)"
    label 'power_analysis'
    publishDir "${params.output_dir}/${params.study_name}/02_Power_analysis", mode: 'copy'
    
    input:
    tuple val(meta), path(original_counts), path(corrected_counts), path(corrected_metadata), path(config_yaml)
    
    output:
    path "*.csv", emit: power_results
    path "power_analysis_log.txt", emit: log
    
    script:
    """
    # Install required R packages if not available
    Rscript -e "
    packages <- c('RNASeqPower', 'DESeq2', 'dplyr', 'yaml', 'argparse')
    new_packages <- packages[!(packages %in% installed.packages()[,'Package'])]
    if(length(new_packages)) {
        if (!requireNamespace('BiocManager', quietly = TRUE)) {
            install.packages('BiocManager', repos='https://cran.r-project.org/')
        }
        BiocManager::install(new_packages, ask = FALSE)
    }
    "
    
    # Run power analysis (CSV outputs only)
    Rscript ${projectDir}/bin/power_analysis.R \\
        --original_counts ${original_counts} \\
        --corrected_counts ${corrected_counts} \\
        --metadata ${corrected_metadata} \\
        --config ${config_yaml} \\
        --output_dir .
    """
    
    stub:
    """
    touch power_analysis_results.csv
    touch gene_specific_power_analysis.csv
    touch qq_analysis_summary.csv
    touch power_analysis_log.txt
    """
}
