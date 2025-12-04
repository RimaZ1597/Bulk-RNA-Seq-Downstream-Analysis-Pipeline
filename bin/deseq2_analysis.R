#!/usr/bin/env Rscript

# DESeq2 Differential Expression Analysis - CSV Output Only (Multi-Dataset Support)
# Usage: Rscript deseq2_analysis.R --corrected_counts <file> --metadata <file> --config <file> --output_dir <dir>

suppressPackageStartupMessages({
    library(DESeq2)       # Differential expression analysis
    library(dplyr)        # Data manipulation
    library(yaml)         # Configuration parsing
    library(argparse)     # Command line arguments
})

# Parse command line arguments
parser <- ArgumentParser(description='DESeq2 differential expression analysis - CSV outputs only with multi-dataset support')
parser$add_argument('--corrected_counts', required=TRUE, help='ComBat-seq corrected count matrix')
parser$add_argument('--metadata', required=TRUE, help='Corrected metadata file')
parser$add_argument('--config', required=TRUE, help='Configuration YAML file')
parser$add_argument('--output_dir', default='.', help='Output directory')
args <- parser$parse_args()

cat("=== DESeq2 Differential Expression Analysis (CSV Only - Multi-Dataset Support) ===\n")
cat("Parameters:\n")
cat("  Corrected counts:", args$corrected_counts, "\n")
cat("  Metadata:", args$metadata, "\n")
cat("  Config:", args$config, "\n")
cat("  Output directory:", args$output_dir, "\n\n")

# Load configuration
cat("Loading configuration...\n")
tryCatch({
    config <- yaml::read_yaml(args$config)
    primary_factor <- config$experimental_design$primary_factor
    additional_factors <- config$experimental_design$additional_factors
    comparisons <- config$comparisons
    analysis_params <- config$analysis_parameters
    
    cat("✓ Configuration loaded successfully\n")
}, error = function(e) {
    stop("Failed to load configuration: ", e$message)
})

# Handle additional factors more robustly
if (is.null(additional_factors) || length(additional_factors) == 0 || 
    (length(additional_factors) == 1 && additional_factors == "")) {
    additional_factors <- NULL
    cat("Additional factors: None\n")
} else {
    # Filter out empty strings
    additional_factors <- additional_factors[additional_factors != ""]
    if (length(additional_factors) == 0) {
        additional_factors <- NULL
    }
    cat("Additional factors:", ifelse(is.null(additional_factors), "None", paste(additional_factors, collapse=", ")), "\n")
}

cat("Primary factor:", primary_factor, "\n")
cat("Total comparisons:", length(comparisons), "\n")

# Load data
cat("Loading corrected data...\n")
tryCatch({
    corrected_counts <- read.csv(args$corrected_counts, row.names = 1, check.names = FALSE)
    sample_df_cleaned <- read.csv(args$metadata, row.names = 1, check.names = FALSE)
    
    corrected_counts <- round(as.matrix(corrected_counts))
    sample_df_cleaned <- sample_df_cleaned[colnames(corrected_counts), ]
    
    cat("✓ Data loaded successfully\n")
    cat("  Count matrix dimensions:", paste(dim(corrected_counts), collapse=" x "), "\n")
    cat("  Metadata dimensions:", paste(dim(sample_df_cleaned), collapse=" x "), "\n")
}, error = function(e) {
    stop("Failed to load data: ", e$message)
})

#####################################
## Helper Functions
#####################################

parse_comparison <- function(comparison_str, sample_df) {
    parts <- strsplit(comparison_str, " vs ")[[1]]
    if (length(parts) != 2) {
        stop("Invalid comparison format: ", comparison_str)
    }
    group1_str <- parts[1]
    group2_str <- parts[2]
    
    parse_group <- function(group_str) {
        if (grepl("Combined_", group_str)) {
            group_str <- gsub("Combined_", "", group_str)
        }
        conditions <- strsplit(group_str, ",")[[1]]
        condition_list <- list()
        for (cond in conditions) {
            if (grepl(":", cond)) {
                cond_parts <- strsplit(cond, ":")[[1]]
                condition_list[[cond_parts[1]]] <- cond_parts[2]
            }
        }
        return(condition_list)
    }
    
    group1_conditions <- parse_group(group1_str)
    group2_conditions <- parse_group(group2_str)
    
    return(list(
        group1 = group1_conditions,
        group2 = group2_conditions,
        name = paste(group1_str, "vs", group2_str)
    ))
}

filter_samples_by_conditions <- function(sample_df, conditions) {
    if (length(conditions) == 0) return(rep(FALSE, nrow(sample_df)))
    filter_logic <- rep(TRUE, nrow(sample_df))
    for (factor_name in names(conditions)) {
        if (factor_name %in% colnames(sample_df)) {
            factor_filter <- sample_df[[factor_name]] == conditions[[factor_name]]
            filter_logic <- filter_logic & factor_filter
        }
    }
    return(filter_logic)
}

#####################################
## Main DESeq2 Analysis
#####################################

cat("\n=== Creating DESeq2 Dataset ===\n")

design_factors <- c(primary_factor)
if (!is.null(additional_factors) && length(additional_factors) > 0) {
    design_factors <- c(design_factors, additional_factors)
}

design_formula <- as.formula(paste("~", paste(design_factors, collapse = " + ")))
cat("Design formula:", deparse(design_formula), "\n")

# Validate that design factors exist in metadata
missing_factors <- design_factors[!design_factors %in% colnames(sample_df_cleaned)]
if (length(missing_factors) > 0) {
    stop("Design factors not found in metadata: ", paste(missing_factors, collapse=", "))
}

# Create DESeq2 dataset
tryCatch({
    dds <- DESeqDataSetFromMatrix(
        countData = corrected_counts,
        colData = sample_df_cleaned,
        design = design_formula
    )
    
    keep <- rowSums(counts(dds)) >= analysis_params$min_counts
    dds <- dds[keep, ]
    
    cat("After filtering: retained", nrow(dds), "genes\n")
}, error = function(e) {
    stop("Failed to create DESeq2 dataset: ", e$message)
})

cat("\n=== Running DESeq2 Analysis ===\n")
tryCatch({
    dds <- DESeq(dds)
    normalized_counts <- counts(dds, normalized = TRUE)
    cat("✓ DESeq2 analysis completed\n")
}, error = function(e) {
    stop("DESeq2 analysis failed: ", e$message)
})

#####################################
## Process Each Comparison Separately
#####################################

cat("\n=== Running Differential Expression Comparisons ===\n")

successful_comparisons <- 0
failed_comparisons <- 0

for (i in 1:length(comparisons)) {
    comparison <- comparisons[i]
    cat("Processing comparison", i, "of", length(comparisons), ":", comparison, "\n")
    
    # Create safe comparison name for folder
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", comparison)
    safe_name <- gsub("__+", "_", safe_name)  # Remove multiple underscores
    safe_name <- gsub("^_|_$", "", safe_name)  # Remove leading/trailing underscores
    
    # Create contrast-specific output directory
    contrast_dir <- file.path(args$output_dir, "03_Differential_Expression_&_Pathway_Analysis", safe_name)
    deseq_output_dir <- file.path(contrast_dir, "Differential_Expression_Analysis")
    dir.create(deseq_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    tryCatch({
        parsed_comp <- parse_comparison(comparison, sample_df_cleaned)
        group1_samples <- filter_samples_by_conditions(sample_df_cleaned, parsed_comp$group1)
        group2_samples <- filter_samples_by_conditions(sample_df_cleaned, parsed_comp$group2)
        
        # Fixed syntax error: added comparison operator
        if (sum(group1_samples) == 0 || sum(group2_samples) == 0) {
            cat("  Warning: No samples found for one or both groups. Skipping.\n")
            failed_comparisons <- failed_comparisons + 1
            next
        }
        
        cat("  Group 1 samples:", sum(group1_samples), "\n")
        cat("  Group 2 samples:", sum(group2_samples), "\n")
        
        comparison_samples <- group1_samples | group2_samples
        dds_subset <- dds[, comparison_samples]
        comparison_factor <- ifelse(group1_samples[comparison_samples], "Group1", "Group2")
        comparison_factor <- factor(comparison_factor, levels = c("Group2", "Group1"))
        colData(dds_subset)$comparison_group <- comparison_factor
        design(dds_subset) <- ~ comparison_group
        dds_subset <- DESeq(dds_subset)
        
        res <- results(dds_subset,
                      contrast = c("comparison_group", "Group1", "Group2"),
                      alpha = analysis_params$padj_threshold)
        
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        
        # Save DESeq2 results (CSV only)
        write.csv(res_df,
                 file.path(deseq_output_dir, "DESeq2_results.csv"),
                 row.names = FALSE)
        
        # Save significant genes
        sig_genes <- res_df[!is.na(res_df$padj) & 
                           res_df$padj < analysis_params$padj_threshold &
                           abs(res_df$log2FoldChange) >= analysis_params$logfc_threshold, ]
        
        write.csv(sig_genes,
                 file.path(deseq_output_dir, "significant_genes.csv"),
                 row.names = FALSE)
        
        # Save normalized counts for this comparison
        comparison_normalized <- normalized_counts[, comparison_samples]
        write.csv(comparison_normalized,
                 file.path(deseq_output_dir, "normalized_counts.csv"))
        
        # Create summary statistics
        upregulated <- sum(sig_genes$log2FoldChange > 0, na.rm = TRUE)
        downregulated <- sum(sig_genes$log2FoldChange < 0, na.rm = TRUE)
        
        summary_stats <- data.frame(
            Comparison = comparison,
            Total_Genes = nrow(res_df),
            Significant_Genes = nrow(sig_genes),
            Upregulated = upregulated,
            Downregulated = downregulated,
            P_value_threshold = analysis_params$pvalue_threshold,
            Adjusted_P_value_threshold = analysis_params$padj_threshold,
            LogFC_threshold = analysis_params$logfc_threshold,
            Group1_Samples = sum(group1_samples),
            Group2_Samples = sum(group2_samples),
            stringsAsFactors = FALSE
        )
        
        write.csv(summary_stats,
                 file.path(deseq_output_dir, "analysis_summary.csv"),
                 row.names = FALSE)
        
        cat("  Results: ", nrow(sig_genes), "significant genes (", upregulated, "up,", downregulated, "down)\n")
        cat("  Saved to:", deseq_output_dir, "\n")
        
        successful_comparisons <- successful_comparisons + 1
        
    }, error = function(e) {
        cat("  Error processing comparison:", e$message, "\n")
        failed_comparisons <- failed_comparisons + 1
    })
}

# Create comprehensive log
log_content <- paste(
    "DESeq2 Analysis Summary Report (CSV Output Only)",
    "===============================================",
    "",
    paste("Timestamp:", Sys.time()),
    paste("Design formula:", deparse(design_formula)),
    paste("Total comparisons requested:", length(comparisons)),
    paste("Successful comparisons:", successful_comparisons),
    paste("Failed comparisons:", failed_comparisons),
    paste("Total genes analyzed:", nrow(dds)),
    paste("Primary factor:", primary_factor),
    paste("Additional factors:", ifelse(is.null(additional_factors), "None", paste(additional_factors, collapse=", "))),
    "",
    "Analysis Parameters:",
    paste("- P-value threshold:", analysis_params$pvalue_threshold),
    paste("- Adjusted p-value threshold:", analysis_params$padj_threshold),
    paste("- Log fold-change threshold:", analysis_params$logfc_threshold),
    paste("- Minimum counts filter:", analysis_params$min_counts),
    "",
    "Output Structure:",
    "03_Differential_Expression_&_Pathway_Analysis/",
    "├── [Contrast_Name]/",
    "│   └── Differential_Expression_Analysis/",
    "│       ├── DESeq2_results.csv",
    "│       ├── significant_genes.csv",
    "│       ├── normalized_counts.csv",
    "│       └── analysis_summary.csv",
    "",
    "Note: All visualizations disabled - CSV outputs only",
    "",
    paste("DESeq2 analysis completed with", successful_comparisons, "successful comparisons!"),
    sep = "\n"
)

writeLines(log_content, file.path(args$output_dir, "deseq2_analysis_log.txt"))

cat("\n=== DESeq2 Analysis Completed! ===\n")
cat("Output format: CSV files only (no visualizations)\n")
cat("Results organized by contrast in separate folders\n")
cat("Successful comparisons:", successful_comparisons, "/", length(comparisons), "\n")

if (failed_comparisons > 0) {
    cat(" failed_comparisons, "comparisons failed - check the log for details\n")
}

cat("Results directory: 03_Differential_Expression_&_Pathway_Analysis/\n")
