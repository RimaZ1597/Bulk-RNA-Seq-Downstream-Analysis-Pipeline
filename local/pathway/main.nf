process PATHWAY_ANALYSIS {
    tag "Pathway enrichment analysis (CSV only, no GO)"
    label 'pathway_analysis'
    publishDir "${params.output_dir}/${params.study_name}", mode: 'copy'
    
    input:
    tuple val(meta), path(deseq2_results_dir), path(config_yaml)
    
    output:
    path "03_Differential_Expression_&_Pathway_Analysis/**/Pathway_Analysis/*.csv", emit: pathway_results
    path "pathway_analysis_log.txt", emit: log
    
    script:
    """
    # Install required R packages if not available
    Rscript -e "
    packages <- c('fgsea', 'msigdbr', 'dplyr', 'yaml', 'argparse')
    new_packages <- packages[!(packages %in% installed.packages()[,'Package'])]
    if(length(new_packages)) {
        if (!requireNamespace('BiocManager', quietly = TRUE)) {
            install.packages('BiocManager', repos='https://cran.r-project.org/')
        }
        BiocManager::install(new_packages, ask = FALSE)
    }
    "
    
    # Run pathway analysis (CSV outputs only, no GO, all pathways)
    Rscript ${projectDir}/bin/pathway_analysis.R \\
        --deseq2_results ${deseq2_results_dir} \\
        --config ${config_yaml} \\
        --output_dir .
    """
    
    stub:
    """
    mkdir -p 03_Differential_Expression_&_Pathway_Analysis/test_contrast/Pathway_Analysis
    touch 03_Differential_Expression_&_Pathway_Analysis/test_contrast/Pathway_Analysis/all_pathway_results.csv
    touch 03_Differential_Expression_&_Pathway_Analysis/test_contrast/Pathway_Analysis/pathway_summary.csv
    touch pathway_analysis_log.txt
    """
}
