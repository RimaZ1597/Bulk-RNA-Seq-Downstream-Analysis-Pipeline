// modules/local/pathway/main.nf

nextflow.enable.dsl = 2

/*
========================================================================================
    FGSEA PATHWAY ANALYSIS MODULE
========================================================================================
    Performs pathway analysis on limma differential expression results
    
    Requirements:
    - Container: Bioconductor image with fgsea, msigdbr packages
    - R script: bin/pathway_analysis.R
    - Input: Limma results (per contrast) + subsystem genes + config
    - Output: CSV files (contrast-specific + combined)
========================================================================================
*/

process PATHWAY_ANALYSIS {
    tag "${project_id}:${contrast_name}"
    label 'process_medium'
    
    publishDir "data/${project_id}/analysis_results/pathway_analysis", mode: 'copy', pattern: "*.csv"
    publishDir "data/${project_id}/analysis_results/limma/${safe_contrast_name}/pathway_analysis", 
               mode: 'copy', pattern: "${safe_contrast_name}/*"
    
    input:
    tuple val(project_id), val(contrast_name), val(safe_contrast_name), 
          path(limma_results), path(config_yaml), path(subsystem_genes)
    path pathway_script
    
    output:
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/ranked_gene_list.csv"), emit: ranked_genes
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/humangem_subsystem_pathways.csv"), emit: humangem_all
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/humangem_subsystem_pathways_significant.csv"), 
          emit: humangem_sig, optional: true
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/msigdb_*_pathways.csv"), emit: msigdb_collections
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/msigdb_all_collections_combined.csv"), 
          emit: msigdb_combined
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/msigdb_all_collections_significant.csv"), 
          emit: msigdb_sig, optional: true
    tuple val(project_id), val(contrast_name), 
          path("${safe_contrast_name}/pathway_analysis_summary.csv"), emit: summary
    path "pathway_analysis_${safe_contrast_name}.log", emit: log
    
    script:
    """
    # Create contrast-specific output directory
    mkdir -p ${safe_contrast_name}
    
    # Run pathway analysis
    Rscript ${pathway_script} \\
        --limma_results ${limma_results} \\
        --contrast_name "${contrast_name}" \\
        --subsystem_genes ${subsystem_genes} \\
        --config ${config_yaml} \\
        --output_dir ${safe_contrast_name} \\
        --species "Homo sapiens" \\
        --fdr_threshold 0.05 \\
        --min_size 15 \\
        --max_size 500 \\
        2>&1 | tee pathway_analysis_${safe_contrast_name}.log
    
    # Validate completion
    if [ \$? -ne 0 ]; then
        echo "ERROR: Pathway analysis failed for ${contrast_name}" >&2
        exit 1
    fi
    
    if [ ! -f ${safe_contrast_name}/pathway_analysis_summary.csv ]; then
        echo "ERROR: pathway_analysis_summary.csv not created" >&2
        exit 1
    fi
    
    echo "✓ Pathway analysis completed successfully for ${contrast_name}"
    """
}

process COMBINE_PATHWAY_RESULTS {
    tag "${project_id}"
    label 'process_low'
    
    publishDir "data/${project_id}/analysis_results/pathway_analysis", mode: 'copy'
    
    input:
    tuple val(project_id), path(all_summaries)
    
    output:
    path "all_contrasts_pathway_summary.csv", emit: combined_summary
    path "pathway_analysis_complete.txt", emit: completion_flag
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Debug: List all files to understand the structure
    cat("=== Working Directory Contents ===\\n")
    all_files <- list.files(recursive = TRUE, full.names = TRUE)
    cat(paste(all_files, collapse = "\\n"), "\\n")
    cat("==================================\\n")
    
    # When files are staged by Nextflow, they might be directly in the work directory
    # Try multiple patterns to find the files
    summary_files <- c()
    
    # Pattern 1: Files in subdirectories (original structure)
    summary_files <- c(summary_files, list.files(pattern = "pathway_analysis_summary.csv", recursive = TRUE, full.names = TRUE))
    
    # Pattern 2: Files directly in working directory (flattened by Nextflow)
    direct_files <- list.files(pattern = "pathway_analysis_summary.csv", full.names = TRUE)
    summary_files <- c(summary_files, direct_files)
    
    # Remove duplicates
    summary_files <- unique(summary_files)
    
    if (length(summary_files) == 0) {
        cat("ERROR: No pathway analysis summary files found\\n")
        cat("Available files:\\n")
        cat(paste(list.files(recursive = TRUE), collapse = "\\n"))
        stop("No pathway analysis summary files found")
    }
    
    cat("Found", length(summary_files), "pathway summary files:\\n")
    cat(paste(summary_files, collapse = "\\n"), "\\n")
    
    # Combine all summaries
    all_summaries <- list()
    
    for (i in seq_along(summary_files)) {
        file <- summary_files[i]
        
        # Extract contrast name from file path or use index
        if (grepl("/", file)) {
            # File is in a subdirectory, extract directory name
            contrast_name <- basename(dirname(file))
        } else {
            # File is in root directory, try to extract from filename or use index
            contrast_name <- paste("contrast", i, sep = "_")
        }
        
        cat("Processing:", file, "-> Contrast:", contrast_name, "\\n")
        
        # Read summary file
        df <- read.csv(file, stringsAsFactors = FALSE)
        df\$contrast <- contrast_name
        
        all_summaries[[contrast_name]] <- df
    }
    
    # Combine all data
    combined <- do.call(rbind, all_summaries)
    
    # Reorder columns
    combined <- combined[c("contrast", "analysis_type", "total_pathways", 
                          "significant_pathways", "fdr_threshold", 
                          "min_pathway_size", "max_pathway_size")]
    
    # Save combined summary
    write.csv(combined, "all_contrasts_pathway_summary.csv", row.names = FALSE, quote = FALSE)
    
    cat("✓ Combined pathway summaries from", length(summary_files), "contrasts\\n")
    cat("✓ Total summary rows:", nrow(combined), "\\n")
    
    # Create completion flag
    writeLines(
        c(paste("Pathway analysis completed for", length(summary_files), "contrasts"),
          paste("Total summaries:", nrow(combined)),
          paste("Timestamp:", Sys.time())),
        "pathway_analysis_complete.txt"
    )
    
    cat("🎉 Pathway analysis combination completed successfully!\\n")
    """
}