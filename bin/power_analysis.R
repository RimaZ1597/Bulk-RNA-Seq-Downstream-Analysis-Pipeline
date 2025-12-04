#!/usr/bin/env Rscript

# Power Analysis - CSV Output Only (Multi-Dataset Support with Combined Analysis)
# Usage: Rscript power_analysis.R --original_counts <file> --corrected_counts <file> --metadata <file> --config <file> --output_dir <dir>

suppressPackageStartupMessages({
    library(RNASeqPower)  # For power calculations
    library(DESeq2)       # For p-value generation
    library(dplyr)        # For data manipulation
    library(yaml)         # For configuration
    library(argparse)     # For command line arguments
})

# Parse command line arguments
parser <- ArgumentParser(description='Power analysis - CSV output only with combined dataset analysis')
parser$add_argument('--original_counts', required=TRUE, help='Original count matrix (before correction)')
parser$add_argument('--corrected_counts', required=TRUE, help='ComBat-seq corrected count matrix')
parser$add_argument('--metadata', required=TRUE, help='Corrected metadata file')
parser$add_argument('--config', required=TRUE, help='Configuration YAML file')
parser$add_argument('--output_dir', default='.', help='Output directory')
args <- parser$parse_args()

cat("=== Power Analysis (CSV Output Only - Enhanced with Combined Analysis) ===\n")
cat("Parameters:\n")
cat("  Original counts:", args$original_counts, "\n")
cat("  Corrected counts:", args$corrected_counts, "\n")
cat("  Metadata:", args$metadata, "\n")
cat("  Config:", args$config, "\n")
cat("  Output directory:", args$output_dir, "\n\n")

# Create output directories
power_dir <- file.path(args$output_dir, "02_Power_analysis")
dir.create(power_dir, showWarnings = FALSE, recursive = TRUE)

# Load configuration
config <- yaml::read_yaml(args$config)
batch_factor <- config$experimental_design$batch_factor
primary_factor <- config$experimental_design$primary_factor
additional_factors <- config$experimental_design$additional_factors

# Handle batch factor and additional factors
skip_batch_analysis <- is.null(batch_factor) || batch_factor == "none" || batch_factor == ""
if (skip_batch_analysis) {
    cat("⚠️ No batch factor specified - performing single-dataset analysis\n")
}

# Handle additional factors
if (is.null(additional_factors) || length(additional_factors) == 0 || 
    (length(additional_factors) == 1 && additional_factors == "")) {
    additional_factors <- NULL
}

# Load data
cat("Loading original and corrected data...\n")
tryCatch({
    original_counts <- read.csv(args$original_counts, row.names = 1, check.names = FALSE)
    corrected_counts <- read.csv(args$corrected_counts, row.names = 1, check.names = FALSE)
    sample_df_cleaned <- read.csv(args$metadata, row.names = 1, check.names = FALSE)
    
    original_counts <- as.matrix(original_counts)
    corrected_counts <- as.matrix(corrected_counts)
    
    cat("✓ Data loaded successfully\n")
}, error = function(e) {
    stop("Failed to load data: ", e$message)
})

# Convert categorical columns to factors based on config
cat("Converting categorical variables to factors...\n")
if (!is.null(config$metadata_info)) {
    for (col_name in names(config$metadata_info)) {
        if (col_name %in% colnames(sample_df_cleaned)) {
            col_info <- config$metadata_info[[col_name]]
            if (!is.null(col_info$type) && col_info$type == "categorical") {
                cat("Converting", col_name, "to factor\n")
                sample_df_cleaned[[col_name]] <- as.factor(sample_df_cleaned[[col_name]])
            }
        }
    }
}

# Also ensure batch factor and other design factors are factors
design_factors <- c(primary_factor)
if (!is.null(additional_factors)) {
    design_factors <- c(design_factors, additional_factors)
}
design_factors <- design_factors[design_factors %in% colnames(sample_df_cleaned)]

for (factor_name in design_factors) {
    if (!is.factor(sample_df_cleaned[[factor_name]])) {
        cat("Converting design factor", factor_name, "to factor\n")
        sample_df_cleaned[[factor_name]] <- as.factor(sample_df_cleaned[[factor_name]])
    }
}

#####################################
## Helper Functions
#####################################

power_analysis_summary <- function(mean_depth, cv, n_samples, label) {
    effect_sizes <- c(1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
    powers <- c()
    for (effect in effect_sizes) {
        power <- RNASeqPower::rnapower(
            depth = mean_depth,
            cv = cv,
            effect = effect,
            alpha = 0.05,
            n = n_samples
        )
        powers <- c(powers, power)
    }
    data.frame(
        Effect_Size = effect_sizes,
        Power = powers,
        Sample_Size = n_samples,
        Mean_Depth = mean_depth,
        CV = cv,
        Label = label,
        stringsAsFactors = FALSE
    )
}

analyze_gene_power_csv <- function(count_data1, count_data2, depth1, depth2, n_samples1, n_samples2, effect_sizes) {
    # Combine data for analysis
    all_data <- cbind(count_data1, count_data2)
    # Calculate CV for each gene
    gene_results <- data.frame(
        gene = rownames(all_data),
        stringsAsFactors = FALSE
    )
    # Calculate CV values
    cv_values <- apply(all_data, 1, function(x) {
        if(all(is.na(x)) || all(x == 0)) return(NA)
        cv <- sd(x)/mean(x)
        if(is.infinite(cv) || is.nan(cv)) return(NA)
        return(cv)
    })
    gene_results$cv <- cv_values
    gene_results$mean_expression <- rowMeans(all_data, na.rm = TRUE)
    gene_results$total_counts <- rowSums(all_data, na.rm = TRUE)
    # Calculate power for each effect size
    for(effect in effect_sizes) {
        power_values <- sapply(cv_values, function(cv) {
            if(is.na(cv)) return(NA)
            tryCatch({
                RNASeqPower::rnapower(
                    depth = mean(c(depth1, depth2)),
                    cv = cv,
                    effect = effect,
                    alpha = 0.05,
                    n = mean(c(n_samples1, n_samples2))
                )
            }, error = function(e) NA)
        })
        gene_results[[paste0("power_", effect, "fold")]] <- power_values
    }
    # Remove genes with NA values
    gene_results_clean <- gene_results[complete.cases(gene_results), ]
    return(gene_results_clean)
}

generate_qq_analysis_csv <- function(count_data, metadata, design_formula, analysis_name) {
    tryCatch({
        cat("Performing DESeq2 analysis for:", analysis_name, "\n")
        # Create DESeq2 dataset
        dds <- DESeqDataSetFromMatrix(
            countData = round(count_data),
            colData = metadata,
            design = design_formula
        )
        # Filter low-count genes
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]
        # Run DESeq2
        dds <- DESeq(dds)
        res <- results(dds)
        # Extract p-values and statistics
        pvals <- res$pvalue[!is.na(res$pvalue)]
        if(length(pvals) > 0) {
            # Calculate QQ plot statistics
            n <- length(pvals)
            expected_pvals <- (1:n) / (n + 1)
            observed_pvals <- sort(pvals)
            qq_data <- data.frame(
                rank = 1:n,
                expected_pval = expected_pvals,
                observed_pval = observed_pvals,
                expected_neg_log10 = -log10(expected_pvals),
                observed_neg_log10 = -log10(observed_pvals),
                analysis = analysis_name,
                stringsAsFactors = FALSE
            )
            # Calculate lambda (genomic inflation factor)
            lambda <- qchisq(median(pvals, na.rm = TRUE), 1, lower.tail = FALSE) / qchisq(0.5, 1)
            return(list(
                qq_data = qq_data,
                summary = data.frame(
                    analysis = analysis_name,
                    total_genes = length(pvals),
                    significant_genes_05 = sum(pvals < 0.05, na.rm = TRUE),
                    significant_genes_01 = sum(pvals < 0.01, na.rm = TRUE),
                    mean_pvalue = mean(pvals, na.rm = TRUE),
                    median_pvalue = median(pvals, na.rm = TRUE),
                    lambda = lambda,
                    stringsAsFactors = FALSE
                )
            ))
        } else {
            return(NULL)
        }
    }, error = function(e) {
        cat("Warning: QQ analysis failed for", analysis_name, ":", e$message, "\n")
        return(NULL)
    })
}

#####################################
## Main Power Analysis
#####################################

cat("\n=== Starting Power Analysis ===\n")

# Initialize variables for combined analysis
all_power_results_list <- list()

# Check if we should do batch-specific analysis
if (!skip_batch_analysis && batch_factor %in% colnames(sample_df_cleaned)) {
    unique_batches <- unique(sample_df_cleaned[[batch_factor]])
    cat("Found batches:", paste(unique_batches, collapse=", "), "\n")
    
    # Split data by batches
    batch_data <- list()
    batch_stats <- list()
    power_results_list <- list()
    
    for (batch in unique_batches) {
        batch_samples <- sample_df_cleaned[[batch_factor]] == batch
        batch_data[[as.character(batch)]] <- corrected_counts[, batch_samples]
        
        # Calculate statistics
        mean_depth <- mean(colSums(batch_data[[as.character(batch)]])) / nrow(batch_data[[as.character(batch)]])
        cv <- mean(apply(batch_data[[as.character(batch)]], 1, function(x) sd(x)/mean(x)), na.rm=TRUE)
        n_samples <- ncol(batch_data[[as.character(batch)]])
        
        batch_stats[[as.character(batch)]] <- list(
            mean_depth = mean_depth,
            cv = cv,
            n_samples = n_samples
        )
        
        cat(paste("Batch", batch, "- Samples:", n_samples, "Mean depth:", format(mean_depth, scientific=TRUE), "CV:", round(cv, 3)), "\n")
        
        # Perform power analysis for each batch
        power_results <- power_analysis_summary(
            batch_stats[[as.character(batch)]]$mean_depth,
            batch_stats[[as.character(batch)]]$cv,
            batch_stats[[as.character(batch)]]$n_samples,
            paste("Batch", batch, "(Corrected)")
        )
        power_results_list[[as.character(batch)]] <- power_results
        
        # Save individual results (CSV only) - fixed path
        write.csv(
            power_results,
            file.path(power_dir, paste0("power_analysis_batch", batch, "_corrected.csv")),
            row.names = FALSE
        )
    }
    
    # Add batch results to combined list
    all_power_results_list <- power_results_list
    
    cat("✓ Batch-specific power analysis completed\n")
} else {
    cat("Performing single-dataset analysis (no batch factor)\n")
}

#####################################
## Combined Dataset Power Analysis (Enhanced)
#####################################

cat("\n=== Combined Dataset Power Analysis (All Plates/Batches Together) ===\n")

# Always perform combined analysis regardless of batch factor
total_samples <- ncol(corrected_counts)
total_genes <- nrow(corrected_counts)
mean_depth_combined <- mean(colSums(corrected_counts)) / nrow(corrected_counts)
cv_combined <- mean(apply(corrected_counts, 1, function(x) {
    if(length(unique(x)) > 1) sd(x)/mean(x) else 0
}), na.rm=TRUE)

cat("Combined Dataset Statistics:\n")
cat("  Total samples (all plates/batches):", total_samples, "\n")
cat("  Total genes:", total_genes, "\n")
cat("  Mean sequencing depth:", format(mean_depth_combined, scientific=TRUE), "\n")
cat("  Overall coefficient of variation:", round(cv_combined, 3), "\n")

# Calculate power across multiple effect sizes for the entire dataset
effect_sizes <- c(1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0)
power_combined_detailed <- power_analysis_summary(
    mean_depth_combined,
    cv_combined,
    total_samples,
    "Combined_All_Plates_Batches"
)

# Add additional statistics
power_combined_detailed$Dataset_Type <- "Combined_All_Samples"
power_combined_detailed$Analysis_Date <- as.character(Sys.time())

# Save combined dataset power analysis
write.csv(
    power_combined_detailed,
    file.path(power_dir, "power_analysis_combined_all_dataset.csv"),
    row.names = FALSE
)

cat("✓ Combined dataset power analysis completed\n")

# Create detailed summary for different sample size scenarios
cat("\n=== Power Analysis for Different Sample Size Scenarios ===\n")

sample_size_scenarios <- c(10, 20, 30, 50, 100, 200)
sample_size_scenarios <- sample_size_scenarios[sample_size_scenarios <= total_samples]

power_scenarios <- data.frame()
for (n in sample_size_scenarios) {
    scenario_power <- power_analysis_summary(
        mean_depth_combined,
        cv_combined,
        n,
        paste("Scenario_N", n, sep="_")
    )
    power_scenarios <- rbind(power_scenarios, scenario_power)
}

# Add current dataset size if not already included
if (!total_samples %in% sample_size_scenarios) {
    current_power <- power_analysis_summary(
        mean_depth_combined,
        cv_combined,
        total_samples,
        paste("Current_Dataset_N", total_samples, sep="_")
    )
    power_scenarios <- rbind(power_scenarios, current_power)
}

write.csv(
    power_scenarios,
    file.path(power_dir, "power_analysis_sample_size_scenarios.csv"),
    row.names = FALSE
)

cat("✓ Sample size scenario analysis completed\n")

#####################################
## Whole Dataset Gene-Level Power Analysis
#####################################

cat("\n=== Whole Dataset Gene-Level Power Analysis ===\n")

# Calculate gene-specific power for the entire combined dataset
gene_power_combined <- data.frame(
    gene = rownames(corrected_counts),
    stringsAsFactors = FALSE
)

# Calculate per-gene statistics across all samples
gene_means <- rowMeans(corrected_counts, na.rm = TRUE)
gene_sds <- apply(corrected_counts, 1, sd, na.rm = TRUE)
gene_cvs <- gene_sds / gene_means
gene_cvs[is.infinite(gene_cvs) | is.nan(gene_cvs)] <- NA

gene_power_combined$mean_expression <- gene_means
gene_power_combined$cv <- gene_cvs
gene_power_combined$total_counts <- rowSums(corrected_counts, na.rm = TRUE)

# Calculate power for each gene at different effect sizes
effect_sizes_gene <- c(1.5, 2.0, 2.5, 3.0)
for (effect in effect_sizes_gene) {
    power_values <- sapply(gene_cvs, function(cv) {
        if (is.na(cv)) return(NA)
        tryCatch({
            RNASeqPower::rnapower(
                depth = mean_depth_combined,
                cv = cv,
                effect = effect,
                alpha = 0.05,
                n = total_samples
            )
        }, error = function(e) NA)
    })
    gene_power_combined[[paste0("power_", effect, "fold")]] <- power_values
}

# Remove genes with missing data
gene_power_combined_clean <- gene_power_combined[complete.cases(gene_power_combined), ]

# Save whole dataset gene-level power analysis
write.csv(
    gene_power_combined_clean,
    file.path(power_dir, "gene_power_whole_dataset.csv"),
    row.names = FALSE
)

# Create summary statistics for whole dataset
whole_dataset_summary <- data.frame(
    analysis_type = "Whole_Dataset_Combined",
    total_genes = nrow(gene_power_combined_clean),
    total_samples = total_samples,
    mean_cv = mean(gene_power_combined_clean$cv, na.rm = TRUE),
    median_cv = median(gene_power_combined_clean$cv, na.rm = TRUE),
    mean_expression = mean(gene_power_combined_clean$mean_expression, na.rm = TRUE),
    genes_power_80_1_5fold = sum(gene_power_combined_clean$power_1.5fold >= 0.8, na.rm = TRUE),
    genes_power_80_2fold = sum(gene_power_combined_clean$power_2fold >= 0.8, na.rm = TRUE),
    genes_power_80_2_5fold = sum(gene_power_combined_clean$power_2.5fold >= 0.8, na.rm = TRUE),
    percent_genes_power_80_2fold = round(100 * sum(gene_power_combined_clean$power_2fold >= 0.8, na.rm = TRUE) / nrow(gene_power_combined_clean), 2),
    stringsAsFactors = FALSE
)

write.csv(
    whole_dataset_summary,
    file.path(power_dir, "whole_dataset_power_summary.csv"),
    row.names = FALSE
)

cat("✓ Whole dataset gene-level power analysis completed\n")
cat("  Genes with adequate power (≥80%) for 2-fold change:",
    whole_dataset_summary$genes_power_80_2fold,
    "(", whole_dataset_summary$percent_genes_power_80_2fold, "%)\n")

#####################################
## Power Comparison: Individual vs Combined
#####################################

if (exists("batch_data") && length(batch_data) > 1) {
    cat("\n=== Power Comparison: Individual Batches vs Combined Dataset ===\n")
    
    # Create comparison table
    power_comparison <- data.frame()
    
    # Add individual batch results
    for (batch_name in names(batch_stats)) {
        batch_power_2fold <- power_analysis_summary(
            batch_stats[[batch_name]]$mean_depth,
            batch_stats[[batch_name]]$cv,
            batch_stats[[batch_name]]$n_samples,
            paste("Individual_Batch", batch_name)
        )
        
        comparison_row <- data.frame(
            analysis_type = paste("Individual_Batch", batch_name),
            sample_size = batch_stats[[batch_name]]$n_samples,
            power_1_5fold = batch_power_2fold$Power[batch_power_2fold$Effect_Size == 1.5],
            power_2fold = batch_power_2fold$Power[batch_power_2fold$Effect_Size == 2.0],
            power_2_5fold = batch_power_2fold$Power[batch_power_2fold$Effect_Size == 2.5],
            mean_depth = batch_stats[[batch_name]]$mean_depth,
            cv = batch_stats[[batch_name]]$cv,
            stringsAsFactors = FALSE
        )
        power_comparison <- rbind(power_comparison, comparison_row)
    }
    
    # Add combined dataset result
    combined_comparison_row <- data.frame(
        analysis_type = "Combined_All_Batches",
        sample_size = total_samples,
        power_1_5fold = power_combined_detailed$Power[power_combined_detailed$Effect_Size == 1.5],
        power_2fold = power_combined_detailed$Power[power_combined_detailed$Effect_Size == 2.0],
        power_2_5fold = power_combined_detailed$Power[power_combined_detailed$Effect_Size == 2.5],
        mean_depth = mean_depth_combined,
        cv = cv_combined,
        stringsAsFactors = FALSE
    )
    power_comparison <- rbind(power_comparison, combined_comparison_row)
    
    write.csv(
        power_comparison,
        file.path(power_dir, "power_comparison_individual_vs_combined.csv"),
        row.names = FALSE
    )
    
    cat("✓ Power comparison analysis completed\n")
    cat("Power improvement by combining all batches:\n")
    cat("  Individual batches (avg):", round(mean(power_comparison$power_2fold[power_comparison$analysis_type != "Combined_All_Batches"]), 3), "\n")
    cat("  Combined dataset:", round(combined_comparison_row$power_2fold, 3), "\n")
}

# Combine all power results - fixed variable references
all_power_results_list[["Combined_Dataset"]] <- power_combined_detailed

if (length(all_power_results_list) > 0) {
    all_power_results <- do.call(rbind, all_power_results_list)
    write.csv(
        all_power_results,
        file.path(power_dir, "all_power_results.csv"),
        row.names = FALSE
    )
}

#####################################
## Gene-Specific Power Analysis (CSV Only) - For Batch Comparison
#####################################

cat("\n=== Gene-Specific Power Analysis (Batch Comparison) ===\n")

if (exists("batch_data") && length(batch_data) >= 2) {
    batch_names <- names(batch_data)
    batch1 <- batch_names[1]
    batch2 <- batch_names[2]
    effect_sizes <- c(1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
    
    gene_analysis <- analyze_gene_power_csv(
        count_data1 = batch_data[[batch1]],
        count_data2 = batch_data[[batch2]],
        depth1 = batch_stats[[batch1]]$mean_depth,
        depth2 = batch_stats[[batch2]]$mean_depth,
        n_samples1 = batch_stats[[batch1]]$n_samples,
        n_samples2 = batch_stats[[batch2]]$n_samples,
        effect_sizes = effect_sizes
    )
    
    # Save gene-specific results (CSV only) - fixed path
    write.csv(
        gene_analysis,
        file.path(power_dir, "gene_specific_power_analysis_batch_comparison.csv"),
        row.names = FALSE
    )
    
    # Create summary statistics
    summary_stats <- data.frame(
        total_genes = nrow(gene_analysis),
        mean_cv = mean(gene_analysis$cv, na.rm = TRUE),
        median_cv = median(gene_analysis$cv, na.rm = TRUE),
        mean_expression = mean(gene_analysis$mean_expression, na.rm = TRUE),
        power_80_percent_1_5fold = sum(gene_analysis$power_1.5fold >= 0.8, na.rm = TRUE),
        power_80_percent_2fold = sum(gene_analysis$power_2fold >= 0.8, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    
    write.csv(
        summary_stats,
        file.path(power_dir, "gene_power_summary_batch_comparison.csv"),
        row.names = FALSE
    )
    
    cat("✓ Gene-specific power analysis (batch comparison) completed\n")
} else {
    cat("Skipping gene-specific batch comparison - need at least 2 batches\n")
}

#####################################
## QQ Plot Analysis (CSV Output Only)
#####################################

cat("\n=== QQ Plot Analysis (CSV Output) ===\n")

# Build design formula more robustly
design_factors <- c(primary_factor)
if (!is.null(additional_factors) && length(additional_factors) > 0) {
    design_factors <- c(design_factors, additional_factors)
}
design_factors <- design_factors[design_factors %in% colnames(sample_df_cleaned)]

if (length(design_factors) > 0) {
    design_formula <- as.formula(paste("~", paste(design_factors, collapse = " + ")))
    cat("Using design formula:", deparse(design_formula), "\n")
    
    # Generate QQ analysis for both datasets
    original_qq <- generate_qq_analysis_csv(
        count_data = original_counts,
        metadata = sample_df_cleaned,
        design_formula = design_formula,
        analysis_name = "Before_Batch_Correction"
    )
    
    corrected_qq <- generate_qq_analysis_csv(
        count_data = corrected_counts,
        metadata = sample_df_cleaned,
        design_formula = design_formula,
        analysis_name = "After_ComBat_Correction"
    )
    
    # Save QQ analysis results - fixed paths
    if (!is.null(original_qq)) {
        write.csv(
            original_qq$qq_data,
            file.path(power_dir, "qq_data_before_correction.csv"),
            row.names = FALSE
        )
        qq_summary <- original_qq$summary
    } else {
        qq_summary <- data.frame(
            analysis = "Before_Batch_Correction",
            total_genes = NA,
            significant_genes_05 = NA,
            significant_genes_01 = NA,
            mean_pvalue = NA,
            median_pvalue = NA,
            lambda = NA,
            stringsAsFactors = FALSE
        )
    }
    
    if (!is.null(corrected_qq)) {
        write.csv(
            corrected_qq$qq_data,
            file.path(power_dir, "qq_data_after_correction.csv"),
            row.names = FALSE
        )
        qq_summary <- rbind(qq_summary, corrected_qq$summary)
    } else {
        qq_summary <- rbind(qq_summary, data.frame(
            analysis = "After_ComBat_Correction",
            total_genes = NA,
            significant_genes_05 = NA,
            significant_genes_01 = NA,
            mean_pvalue = NA,
            median_pvalue = NA,
            lambda = NA,
            stringsAsFactors = FALSE
        ))
    }
    
    # Save QQ summary
    write.csv(
        qq_summary,
        file.path(power_dir, "qq_analysis_summary.csv"),
        row.names = FALSE
    )
    
    cat("✓ QQ plot analysis completed\n")
} else {
    cat("⚠️ No valid design factors found - skipping QQ analysis\n")
}

#####################################
## Summary Report
#####################################

cat("\n=== Creating Summary Report ===\n")

# Create comprehensive log
log_content <- paste(
    "Power Analysis Summary Report (CSV Output Only - Enhanced)",
    "=========================================================",
    "",
    paste("Timestamp:", Sys.time()),
    paste("Total genes analyzed:", nrow(corrected_counts)),
    paste("Total samples:", ncol(corrected_counts)),
    paste("Batches analyzed:", if(exists("unique_batches")) length(unique_batches) else 1),
    paste("Batch factor:", ifelse(skip_batch_analysis, "None", batch_factor)),
    "",
    "Power Analysis Results:",
    if(exists("batch_stats")) {
        paste("Batch-specific power analysis completed for", length(batch_stats), "batches")
    } else {
        "Single dataset power analysis completed"
    },
    paste("Combined dataset power analysis completed with", total_samples, "samples"),
    paste("Gene-level power analysis completed for", nrow(gene_power_combined_clean), "genes"),
    "",
    "Enhanced Features:",
    "- Combined dataset power analysis (all samples together)",
    "- Sample size scenario analysis",
    "- Whole dataset gene-level power calculation",
    "- Individual vs combined power comparison",
    "",
    "QQ Plot Analysis Results:",
    if(exists("original_qq") && !is.null(original_qq)) {
        paste("Before correction - Genes:", original_qq$summary$total_genes,
              "Significant (p<0.05):", original_qq$summary$significant_genes_05,
              "Lambda:", round(original_qq$summary$lambda, 3))
    } else {
        "Before correction - Analysis failed or skipped"
    },
    if(exists("corrected_qq") && !is.null(corrected_qq)) {
        paste("After correction - Genes:", corrected_qq$summary$total_genes,
              "Significant (p<0.05):", corrected_qq$summary$significant_genes_05,
              "Lambda:", round(corrected_qq$summary$lambda, 3))
    } else {
        "After correction - Analysis failed or skipped"
    },
    "",
    "Output Files Generated:",
    "- power_analysis_*.csv: Power analysis results by batch",
    "- power_analysis_combined_all_dataset.csv: Combined dataset power",
    "- power_analysis_sample_size_scenarios.csv: Different sample sizes",
    "- gene_power_whole_dataset.csv: Gene-level power for whole dataset",
    "- whole_dataset_power_summary.csv: Summary statistics",
    "- power_comparison_individual_vs_combined.csv: Comparison analysis",
    "- gene_specific_power_analysis_*.csv: Gene-level power calculations",
    "- qq_data_*.csv: QQ plot data (before/after correction)",
    "- qq_analysis_summary.csv: QQ plot summary statistics",
    "",
    "Note: All visualizations disabled - CSV outputs only",
    "",
    "Enhanced power analysis completed successfully!",
    sep = "\n"
)

writeLines(log_content, file.path(power_dir, "power_analysis_log.txt"))

cat("\n=== Power Analysis Completed Successfully ===\n")
