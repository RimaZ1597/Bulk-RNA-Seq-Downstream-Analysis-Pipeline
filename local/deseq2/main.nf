process DESEQ2_ANALYSIS {
    tag "DESeq2 differential expression (CSV only)"
    label 'deseq2'
    publishDir "${params.output_dir}/${params.study_name}", mode: 'copy'
    
    input:
    tuple val(meta), path(corrected_counts), path(corrected_metadata), path(config_yaml)
    
    output:
    tuple val(meta), path("03_Differential_Expression_&_Pathway_Analysis/"), path(config_yaml), emit: deseq2_results
    path "deseq2_analysis_log.txt", emit: log
    
    script:
    """
    # Install required R packages if not available
    Rscript -e "
    packages <- c('DESeq2', 'dplyr', 'yaml', 'argparse')
    new_packages <- packages[!(packages %in% installed.packages()[,'Package'])]
    if(length(new_packages)) {
        if (!requireNamespace('BiocManager', quietly = TRUE)) {
            install.packages('BiocManager', repos='https://cran.r-project.org/')
        }
        BiocManager::install(new_packages, ask = FALSE)
    }
    "
    
    # Run DESeq2 analysis (CSV outputs only, organized by contrast)
    Rscript ${projectDir}/bin/deseq2_analysis.R \\
        --corrected_counts ${corrected_counts} \\
        --metadata ${corrected_metadata} \\
        --config ${config_yaml} \\
        --output_dir .
    """
    
    stub:
    """
    mkdir -p 03_Differential_Expression_&_Pathway_Analysis/test_contrast/Differential_Expression_Analysis
    touch 03_Differential_Expression_&_Pathway_Analysis/test_contrast/Differential_Expression_Analysis/DESeq2_results.csv
    touch 03_Differential_Expression_&_Pathway_Analysis/test_contrast/Differential_Expression_Analysis/significant_genes.csv
    touch deseq2_analysis_log.txt
    """
}
