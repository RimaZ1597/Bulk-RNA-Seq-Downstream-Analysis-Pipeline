process COMBAT_SEQ {
    tag "ComBat-seq batch correction"
    label 'combat_seq'
    publishDir "${params.output_dir}/${params.study_name}/01_Batch_correction", mode: 'copy'
    
    input:
    tuple val(meta), path(counts_matrix), path(metadata), path(config_yaml)
    
    output:
    tuple val(meta), path("corrected_counts.csv"), path(metadata), emit: corrected_data
    path "batch_correction_log.txt", emit: log
    
    script:
    """
    # Install required R packages if not available
    Rscript -e "
    packages <- c('sva', 'DESeq2', 'yaml', 'argparse')
    new_packages <- packages[!(packages %in% installed.packages()[,'Package'])]
    if(length(new_packages)) {
        if (!requireNamespace('BiocManager', quietly = TRUE)) {
            install.packages('BiocManager', repos='https://cran.r-project.org/')
        }
        BiocManager::install(new_packages, ask = FALSE)
    }
    "
    
    # Run ComBat-seq batch correction (CSV outputs only)
    Rscript ${projectDir}/bin/combat_seq.R \\
        --counts ${counts_matrix} \\
        --metadata ${metadata} \\
        --config ${config_yaml} \\
        --output_dir .
    """
    
    stub:
    """
    touch corrected_counts.csv
    touch batch_correction_log.txt
    """
}
