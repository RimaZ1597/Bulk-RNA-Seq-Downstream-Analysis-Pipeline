#!/usr/bin/env Rscript

###############################################################################
# RNA-seq Power Analysis using RNASeqPower
#
# Based on:
# - Hart et al. "Calculating Sample Size Estimates for RNA Sequencing Data"
# - RNASeqPower Bioconductor package
#
# Theory from package documentation:
#   var(log(y)) ≈ 1/μ + ψ²
#   where μ = depth, ψ = CV
#   
#   Power calculations require:
#   - Sufficient depth (μ ≥ 10) for reliable CV estimation
#   - Filter genes where noise (1/μ) dominates signal (ψ²)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(RNASeqPower)
})

# -----------------------------------------------------------------------------
# Command Line Interface - Base R Implementation (no optparse dependency)
# -----------------------------------------------------------------------------

# Helper function to parse command line arguments
parse_args <- function(args) {
  parsed <- list()
  i <- 1
  
  while (i <= length(args)) {
    arg <- args[i]
    
    if (startsWith(arg, "--")) {
      key <- gsub("^--", "", arg)
      
      if (i < length(args) && !startsWith(args[i + 1], "--")) {
        parsed[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        parsed[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  return(parsed)
}

# Helper function for default values (define before use)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(args)

# Define required parameters
required_params <- c("original_counts", "corrected_counts", "metadata")

# Set parameter defaults
opt$experiment_id <- opt$experiment_id %||% "PROJECT"
opt$outdir <- opt$outdir %||% "."
opt$condition_column <- opt$condition_column %||% "condition"
# batch_column can be NULL
opt$min_count <- as.numeric(opt$min_count %||% 10)
opt$min_prop <- as.numeric(opt$min_prop %||% 0.5)
opt$effect_sizes <- opt$effect_sizes %||% "1.25,1.5,1.75,2"
opt$power_levels <- opt$power_levels %||% "0.8,0.9"
opt$alpha <- as.numeric(opt$alpha %||% 0.05)

# Check that required parameters are provided
missing_params <- required_params[!required_params %in% names(opt)]
if (length(missing_params) > 0) {
  cat("Error: Missing required parameters:", paste(missing_params, collapse = ", "), "\n")
  cat("Usage: Rscript power_analysis.R --original_counts file.csv --corrected_counts file.csv --metadata file.csv [options]\n")
  quit(status = 1)
}

# Parse comma-separated numeric parameters
opt$effect_sizes <- as.numeric(strsplit(opt$effect_sizes, ",")[[1]])
opt$power_levels <- as.numeric(strsplit(opt$power_levels, ",")[[1]])

# -----------------------------------------------------------------------------
# Input Validation
# -----------------------------------------------------------------------------

validate_inputs <- function(opt) {
  errors <- c()
  
  # Check file existence
  if (!file.exists(opt$original_counts)) {
    errors <- c(errors, sprintf("Original counts not found: %s", opt$original_counts))
  }
  if (!file.exists(opt$corrected_counts)) {
    errors <- c(errors, sprintf("Corrected counts not found: %s", opt$corrected_counts))
  }
  if (!file.exists(opt$metadata)) {
    errors <- c(errors, sprintf("Metadata not found: %s", opt$metadata))
  }
  
  # Check parameter validity
  if (opt$min_count < 1) {
    errors <- c(errors, "min_count must be >= 1")
  }
  if (opt$min_prop <= 0 || opt$min_prop > 1) {
    errors <- c(errors, "min_prop must be between 0 and 1")
  }
  if (opt$alpha <= 0 || opt$alpha >= 1) {
    errors <- c(errors, "alpha must be between 0 and 1")
  }
  if (any(opt$effect_sizes <= 1)) {
    errors <- c(errors, "All effect_sizes must be > 1")
  }
  if (any(opt$power_levels <= 0) || any(opt$power_levels >= 1)) {
    errors <- c(errors, "All power_levels must be between 0 and 1")
  }
  
  if (length(errors) > 0) {
    cat("❌ Input validation failed:\n")
    for (err in errors) cat("  - ", err, "\n")
    quit(status = 1)
  }
  
  invisible(TRUE)
}

validate_inputs(opt)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Calculate power statistics for a dataset
#'
#' Implements filtering and CV calculation following RNASeqPower theory:
#' - Filters genes with insufficient counts (noise dominates)
#' - Calculates CV per condition
#' - Returns statistics for power calculations
#'
#' @param count_data Numeric matrix (genes x samples)
#' @param metadata Data frame with sample metadata
#' @param group_name Character string for this analysis group
#' @param condition_col Name of condition column in metadata
#' @param min_count Minimum count threshold
#' @param min_prop Minimum proportion of samples passing threshold
#' @return Data frame with power statistics per condition
calculate_power_stats <- function(count_data, metadata, group_name,
                                   condition_col, min_count, min_prop) {
  
  # Validate condition column exists
  if (!condition_col %in% colnames(metadata)) {
    stop(sprintf("❌ Condition column '%s' not found in metadata", condition_col))
  }
  
  # Align samples between counts and metadata
  common_samples <- intersect(colnames(count_data), rownames(metadata))
  if (length(common_samples) == 0) {
    stop("❌ No overlapping samples between count matrix and metadata")
  }
  
  count_data <- count_data[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]
  condition <- metadata[[condition_col]]
  
  # Validate we have at least 2 conditions
  unique_conditions <- unique(condition)
  if (length(unique_conditions) < 2) {
    stop(sprintf("❌ Need at least 2 conditions. Found: %s", 
                 paste(unique_conditions, collapse = ", ")))
  }
  
  # Filter low-expression genes
  # Theory: Keep genes where μ ≥ min_count in sufficient samples
  # This ensures 1/μ term in var(log(y)) ≈ 1/μ + ψ² is small
  keep <- rowSums(count_data >= min_count) >= ncol(count_data) * min_prop
  filtered_data <- count_data[keep, , drop = FALSE]
  
  n_genes_kept <- sum(keep)
  n_genes_total <- nrow(count_data)
  pct_kept <- 100 * mean(keep)
  
  cat(sprintf("  Filtered: %d/%d genes retained (%.1f%%)\n",
              n_genes_kept, n_genes_total, pct_kept))
  
  if (n_genes_kept < 100) {
    warning(sprintf("⚠️  Only %d genes passed filter. Results may be unreliable.", n_genes_kept))
  }
  
  # First, identify conditions with sufficient samples for power analysis
  sample_counts_per_condition <- table(condition)
  valid_conditions <- names(sample_counts_per_condition[sample_counts_per_condition >= 2])
  insufficient_conditions <- names(sample_counts_per_condition[sample_counts_per_condition < 2])
  
  # Report filtering information
  if (length(insufficient_conditions) > 0) {
    cat("⚠️  Excluding conditions with insufficient samples for power analysis:\n")
    for (cond in insufficient_conditions) {
      cat(sprintf("     • Condition '%s': n = %d (need ≥2)\n", cond, sample_counts_per_condition[cond]))
    }
    cat(sprintf("✓ Valid conditions for analysis: %s (n ≥ 2)\n\n", paste(valid_conditions, collapse = ", ")))
  }
  
  # Calculate CV per condition (only for valid conditions)
  # CV = σ/μ, measures within-group variability
  cv_by_condition <- sapply(valid_conditions, function(cond) {
    samples <- rownames(metadata)[condition == cond]
    
    # Calculate gene-wise CV for this condition
    gene_cvs <- apply(filtered_data[, samples, drop = FALSE], 1, function(x) {
      mu <- mean(x)
      if (mu == 0) return(NA_real_)
      sigma <- sd(x)
      sigma / mu
    })
    
    # Return median CV (robust to outliers)
    median(gene_cvs, na.rm = TRUE)
  })
  
  # Create data frame for valid conditions only
  valid_condition_names <- names(cv_by_condition)
  n_samples_per_valid_condition <- as.numeric(sample_counts_per_condition[valid_condition_names])
  
  # Build results for valid conditions only
  results_list <- list()
  for (i in seq_along(valid_condition_names)) {
    cond_name <- valid_condition_names[i]
    cond_cv <- cv_by_condition[i]
    cond_n <- n_samples_per_valid_condition[i]
    
    results_list[[i]] <- data.frame(
      Group = group_name,
      Condition = cond_name,
      N_Samples = cond_n,
      CV = cond_cv,
      Mean_Depth = mean(colSums(filtered_data)) / nrow(filtered_data),
      Genes_Analyzed = n_genes_kept,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all results
  if (length(results_list) > 0) {
    power_stats <- do.call(rbind, results_list)
  } else {
    warning("⚠️  No conditions found for power analysis")
    power_stats <- data.frame(
      Group = character(0),
      Condition = character(0),
      N_Samples = numeric(0),
      CV = numeric(0),
      Mean_Depth = numeric(0),
      Genes_Analyzed = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(power_stats)
}

#' Calculate actual power achieved with current sample sizes
#'
#' Uses RNASeqPower::rnapower() to calculate power for each
#' combination of condition and effect size
#'
#' @param power_stats Data frame from calculate_power_stats()
#' @param effect_sizes Numeric vector of fold changes
#' @param alpha Significance level
#' @return Data frame with power calculations
calculate_actual_power <- function(power_stats, effect_sizes, alpha = 0.05) {
  
  # Create grid of all combinations
  results <- expand.grid(
    Group = unique(power_stats$Group),
    Condition = power_stats$Condition,
    Effect_Size = effect_sizes,
    stringsAsFactors = FALSE
  )
  
  # Merge with power statistics
  results <- merge(
    results, 
    power_stats[, c("Group", "Condition", "N_Samples", "CV", "Mean_Depth")],
    by = c("Group", "Condition")
  )
  
  # Calculate power for each row using RNASeqPower
  results$Actual_Power <- mapply(
    function(depth, n, cv, effect) {
      # Handle missing values
      if (is.na(cv) || is.na(n) || n < 2) return(NA_real_)
      
      # Call RNASeqPower::rnapower()
      tryCatch({
        rnapower(
          depth = depth,
          n = n,
          cv = cv,
          effect = effect,
          alpha = alpha
        )
      }, error = function(e) {
        warning(sprintf("rnapower failed: depth=%.1f, n=%d, cv=%.3f, effect=%.2f", 
                        depth, n, cv, effect))
        NA_real_
      })
    }, 
    results$Mean_Depth, 
    results$N_Samples, 
    results$CV, 
    results$Effect_Size
  )
  
  return(results)
}

#' Calculate required sample sizes for future studies
#'
#' Uses RNASeqPower::rnapower() to solve for n given target power
#'
#' @param power_stats Data frame from calculate_power_stats()
#' @param effect_sizes Numeric vector of fold changes
#' @param power_levels Numeric vector of target power levels
#' @param alpha Significance level
#' @return Data frame with sample size recommendations
calculate_required_sample_size <- function(power_stats, effect_sizes, 
                                           power_levels, alpha = 0.05) {
  
  # Use median CV and mean depth across all conditions
  median_cv <- median(power_stats$CV, na.rm = TRUE)
  mean_depth <- mean(power_stats$Mean_Depth, na.rm = TRUE)
  
  # Create grid of effect sizes and power levels
  results <- expand.grid(
    Effect_Size = effect_sizes,
    Target_Power = power_levels,
    stringsAsFactors = FALSE
  )
  
  # Calculate required sample size for each combination
  results$Required_N <- mapply(
    function(effect, power_target) {
      tryCatch({
        # Solve for n given power
        n_required <- rnapower(
          depth = mean_depth,
          cv = median_cv,
          effect = effect,
          alpha = alpha,
          power = power_target
        )
        ceiling(n_required)
      }, error = function(e) {
        warning(sprintf("Sample size calculation failed for effect=%.2f, power=%.2f", 
                        effect, power_target))
        NA_integer_
      })
    }, 
    results$Effect_Size, 
    results$Target_Power
  )
  
  results$Median_CV <- median_cv
  results$Mean_Depth <- mean_depth
  
  return(results)
}

# -----------------------------------------------------------------------------
# Main Analysis Pipeline
# -----------------------------------------------------------------------------

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("RNA-seq Power Analysis using RNASeqPower\n")
cat(rep("=", 80), "\n", sep = "")
cat("Experiment ID     :", opt$experiment_id, "\n")
cat("Condition column  :", opt$condition_column, "\n")
cat("Batch column      :", ifelse(is.null(opt$batch_column), "None (global analysis)", opt$batch_column), "\n")
cat("Min count filter  :", opt$min_count, "reads in", opt$min_prop * 100, "% of samples\n")
cat("Effect sizes      :", paste(opt$effect_sizes, collapse = ", "), "\n")
cat("Power levels      :", paste(opt$power_levels, collapse = ", "), "\n")
cat("Alpha level       :", opt$alpha, "\n")
cat(rep("=", 80), "\n\n", sep = "")

# Create output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Read Input Data
# -----------------------------------------------------------------------------

cat("📖 Reading input files...\n")

orig_counts <- read.csv(opt$original_counts, row.names = 1, check.names = FALSE)
corr_counts <- read.csv(opt$corrected_counts, row.names = 1, check.names = FALSE)
metadata <- read.csv(opt$metadata, row.names = 1, check.names = FALSE)

cat(sprintf("  ✓ Original counts  : %d genes × %d samples\n", nrow(orig_counts), ncol(orig_counts)))
cat(sprintf("  ✓ Corrected counts : %d genes × %d samples\n", nrow(corr_counts), ncol(corr_counts)))
cat(sprintf("  ✓ Metadata         : %d samples × %d variables\n", nrow(metadata), ncol(metadata)))

# Align samples across all inputs
common_samples <- Reduce(intersect, list(
  colnames(orig_counts),
  colnames(corr_counts),
  rownames(metadata)
))

if (length(common_samples) == 0) {
  stop("❌ No overlapping samples between counts and metadata")
}

cat(sprintf("  ✓ Common samples   : %d\n\n", length(common_samples)))

# Subset to common samples
orig_counts <- orig_counts[, common_samples, drop = FALSE]
corr_counts <- corr_counts[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# -----------------------------------------------------------------------------
# Step 2: Validate Metadata and Determine Analysis Strategy
# -----------------------------------------------------------------------------

cat("🔬 Analyzing metadata structure...\n")

# Validate condition column
if (!opt$condition_column %in% colnames(metadata)) {
  stop(sprintf(
    "❌ Condition column '%s' not found.\n   Available columns: %s",
    opt$condition_column,
    paste(colnames(metadata), collapse = ", ")
  ))
}

metadata[[opt$condition_column]] <- factor(metadata[[opt$condition_column]])
conditions <- levels(metadata[[opt$condition_column]])

cat(sprintf("  ✓ Conditions found : %s\n", paste(conditions, collapse = ", ")))

# Determine if batch stratification is possible
if (!is.null(opt$batch_column) && opt$batch_column %in% colnames(metadata)) {
  groups <- unique(metadata[[opt$batch_column]])
  strat_var <- opt$batch_column
  cat(sprintf("  ✓ Batch stratification : %d groups (%s)\n\n", 
              length(groups), paste(groups, collapse = ", ")))
} else {
  if (!is.null(opt$batch_column)) {
    cat(sprintf("  ⚠️  Batch column '%s' not found in metadata\n", opt$batch_column))
  }
  cat("  ℹ️  Performing global analysis (no batch stratification)\n\n")
  metadata$analysis_group <- "All_Samples"
  groups <- "All_Samples"
  strat_var <- "analysis_group"
}

# -----------------------------------------------------------------------------
# Step 3: Power Analysis on ORIGINAL Counts
# -----------------------------------------------------------------------------

cat(rep("=", 80), "\n", sep = "")
cat("STEP 1: Power Analysis on ORIGINAL Counts (before batch correction)\n")
cat(rep("=", 80), "\n\n", sep = "")

all_power_stats_original <- list()
all_actual_power_original <- list()

for (grp in groups) {
  cat(rep("─", 80), "\n", sep = "")
  cat(sprintf("📊 Group: %s (ORIGINAL data)\n", grp))
  cat(rep("─", 80), "\n", sep = "")
  
  # Subset data for this group
  group_samples <- rownames(metadata)[metadata[[strat_var]] == grp]
  group_meta <- metadata[group_samples, , drop = FALSE]
  group_counts <- orig_counts[, group_samples, drop = FALSE]
  
  cat(sprintf("  Samples: %d\n", length(group_samples)))
  condition_counts <- table(group_meta[[opt$condition_column]])
  valid_for_power <- condition_counts[condition_counts >= 2]
  invalid_for_power <- condition_counts[condition_counts < 2]
  
  cat(sprintf("  Conditions: %s\n", 
              paste(condition_counts, collapse = " vs ")))
  
  if (length(invalid_for_power) > 0) {
    cat(sprintf("  ⚠️  Note: %d condition(s) excluded from power analysis (n < 2): %s\n",
                length(invalid_for_power), paste(names(invalid_for_power), collapse = ", ")))
    cat(sprintf("  ✓ Analyzing %d valid condition(s): %s\n",
                length(valid_for_power), paste(names(valid_for_power), collapse = ", ")))
  }
  
  # Calculate power statistics
  power_stats <- calculate_power_stats(
    count_data = group_counts,
    metadata = group_meta,
    group_name = grp,
    condition_col = opt$condition_column,
    min_count = opt$min_count,
    min_prop = opt$min_prop
  )
  
  # ✅ SAFETY CHECK: Skip if no valid conditions
  if (nrow(power_stats) == 0) {
    cat("  ⚠️  Skipping power analysis: No conditions with sufficient samples (≥2)\n\n")
    # Store empty results
    all_power_stats_original[[grp]] <- power_stats
    all_actual_power_original[[grp]] <- data.frame(
      Group = character(0),
      stringsAsFactors = FALSE
    )
    next
  }
  
  all_power_stats_original[[grp]] <- power_stats
  
  # Calculate actual power
  actual_power <- calculate_actual_power(
    power_stats = power_stats,
    effect_sizes = opt$effect_sizes,
    alpha = opt$alpha
  )
  actual_power$Dataset <- "Original"
  
  all_actual_power_original[[grp]] <- actual_power
  
  # Display summary
  cat("\n  Power Summary:\n")
  summary_table <- actual_power %>%
    select(Condition, Effect_Size, Actual_Power) %>%
    arrange(Condition, Effect_Size)
  print(summary_table, row.names = FALSE)
  cat("\n")
}

# -----------------------------------------------------------------------------
# Step 4: Power Analysis on CORRECTED Counts
# -----------------------------------------------------------------------------

cat(rep("=", 80), "\n", sep = "")
cat("STEP 2: Power Analysis on CORRECTED Counts (after batch correction)\n")
cat(rep("=", 80), "\n\n", sep = "")

all_power_stats_corrected <- list()
all_actual_power_corrected <- list()

for (grp in groups) {
  cat(rep("─", 80), "\n", sep = "")
  cat(sprintf("📊 Group: %s (CORRECTED data)\n", grp))
  cat(rep("─", 80), "\n", sep = "")
  
  # Subset data for this group
  group_samples <- rownames(metadata)[metadata[[strat_var]] == grp]
  group_meta <- metadata[group_samples, , drop = FALSE]
  group_counts <- corr_counts[, group_samples, drop = FALSE]
  
  cat(sprintf("  Samples: %d\n", length(group_samples)))
  
  # Calculate power statistics
  power_stats <- calculate_power_stats(
    count_data = group_counts,
    metadata = group_meta,
    group_name = grp,
    condition_col = opt$condition_column,
    min_count = opt$min_count,
    min_prop = opt$min_prop
  )
  
  # ✅ SAFETY CHECK: Skip if no valid conditions
  if (nrow(power_stats) == 0) {
    cat("  ⚠️  Skipping power analysis: No conditions with sufficient samples (≥2)\n\n")
    # Store empty results
    all_power_stats_corrected[[grp]] <- power_stats
    all_actual_power_corrected[[grp]] <- data.frame(
      Group = character(0),
      stringsAsFactors = FALSE
    )
    next
  }
  
  all_power_stats_corrected[[grp]] <- power_stats
  
  # Calculate actual power
  actual_power <- calculate_actual_power(
    power_stats = power_stats,
    effect_sizes = opt$effect_sizes,
    alpha = opt$alpha
  )
  actual_power$Dataset <- "Corrected"
  
  all_actual_power_corrected[[grp]] <- actual_power
  
  # Display summary
  cat("\n  Power Summary:\n")
  summary_table <- actual_power %>%
    select(Condition, Effect_Size, Actual_Power) %>%
    arrange(Condition, Effect_Size)
  print(summary_table, row.names = FALSE)
  cat("\n")
}

# -----------------------------------------------------------------------------
# Step 5: Calculate Improvement from Batch Correction
# -----------------------------------------------------------------------------

cat(rep("=", 80), "\n", sep = "")
cat("STEP 3: Batch Correction Impact Analysis\n")
cat(rep("=", 80), "\n\n", sep = "")

# Combine results
power_stats_original <- bind_rows(all_power_stats_original)
power_stats_corrected <- bind_rows(all_power_stats_corrected)
actual_power_original <- bind_rows(all_actual_power_original)
actual_power_corrected <- bind_rows(all_actual_power_corrected)

# ✅ ADDITIONAL SAFETY: Check if combined results are empty
if (nrow(power_stats_corrected) == 0 || nrow(actual_power_corrected) == 0) {
  cat(rep("=", 80), "\n", sep = "")
  cat("⚠️  POWER ANALYSIS CANNOT PROCEED\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("\nReason: No conditions have sufficient samples (≥2) for CV calculation\n\n")
  
  cat("Dataset Summary:\n")
  cat("  Total samples:", length(common_samples), "\n")
  cat("  Conditions:", paste(conditions, collapse = ", "), "\n")
  cat("  Samples per condition:\n")
  print(table(metadata[[opt$condition_column]]))
  cat("\n")
  
  cat("Recommendation:\n")
  cat("  • Each condition needs at least 2 biological replicates\n")
  cat("  • Power analysis requires within-group variance estimation (CV = σ/μ)\n")
  cat("  • With n=1, variance cannot be calculated (division by zero in n-1)\n")
  cat("  • Consider:\n")
  cat("    - Merging similar conditions to increase sample size\n")
  cat("    - Collecting additional biological replicates\n")
  cat("    - Focusing on conditions with adequate replication\n\n")
  
  # Create minimal output files to satisfy Nextflow expectations
  empty_message <- data.frame(
    Message = "Insufficient samples for power analysis. Each condition requires ≥2 samples.",
    Total_Samples = length(common_samples),
    Conditions = paste(conditions, collapse = "; "),
    Samples_Per_Condition = paste(names(table(metadata[[opt$condition_column]])), 
                                   table(metadata[[opt$condition_column]]), 
                                   sep = "=", collapse = "; ")
  )
  
  write.csv(empty_message, file.path(opt$outdir, "power_distribution_original.csv"), row.names = FALSE)
  write.csv(empty_message, file.path(opt$outdir, "power_distribution_corrected.csv"), row.names = FALSE)
  write.csv(empty_message, file.path(opt$outdir, "power_comparison_original_vs_corrected.csv"), row.names = FALSE)
  write.csv(empty_message, file.path(opt$outdir, "power_comparison_combined.csv"), row.names = FALSE)
  write.csv(empty_message, file.path(opt$outdir, "recommended_sample_sizes.csv"), row.names = FALSE)
  write.csv(empty_message, file.path(opt$outdir, "power_comparison_summary.csv"), row.names = FALSE)
  
  # Create log file
  sink(file.path(opt$outdir, "power_analysis_log.txt"))
  cat("Power Analysis Failed - Insufficient Samples\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Experiment ID:", opt$experiment_id, "\n\n")
  cat("Reason: No conditions with ≥2 samples for CV calculation\n\n")
  cat("RNASeqPower Theory:\n")
  cat("  CV = σ/μ where σ = standard deviation, μ = mean\n")
  cat("  σ requires n ≥ 2 (degrees of freedom = n-1)\n")
  cat("  With n=1, CV is mathematically undefined\n\n")
  cat("Sample Distribution:\n")
  print(table(metadata[[opt$condition_column]]))
  cat("\n")
  cat("Recommendation: Collect at least 2 replicates per condition\n")
  sink()
  
  cat("\n⚠️  Created placeholder output files (analysis could not proceed)\n")
  cat("Review power_analysis_log.txt for details\n\n")
  
  # Exit successfully (files created, pipeline can continue)
  quit(status = 0)
}

# Calculate improvement (only if we have valid data)
power_improvement <- actual_power_corrected %>%
  select(Group, Condition, Effect_Size, Actual_Power, N_Samples, CV, Mean_Depth) %>%
  rename(Power_Corrected = Actual_Power, 
         CV_Corrected = CV) %>%
  left_join(
    actual_power_original %>%
      select(Group, Condition, Effect_Size, Actual_Power, CV) %>%
      rename(Power_Original = Actual_Power,
             CV_Original = CV),
    by = c("Group", "Condition", "Effect_Size")
  ) %>%
  mutate(
    Power_Gain = Power_Corrected - Power_Original,
    Power_Gain_Percent = 100 * (Power_Corrected - Power_Original) / pmax(Power_Original, 0.01),
    CV_Reduction = CV_Original - CV_Corrected,
    CV_Reduction_Percent = 100 * (CV_Original - CV_Corrected) / CV_Original
  )

cat("📈 Improvement Summary:\n\n")

improvement_summary <- power_improvement %>%
  group_by(Effect_Size) %>%
  summarise(
    Mean_CV_Reduction_Pct = mean(CV_Reduction_Percent, na.rm = TRUE),
    Mean_Power_Gain = mean(Power_Gain, na.rm = TRUE),
    Mean_Power_Gain_Pct = mean(Power_Gain_Percent, na.rm = TRUE),
    .groups = "drop"
  )

print(improvement_summary, row.names = FALSE)
cat("\n")

# Overall statistics
cat("Overall Impact:\n")
cat(sprintf("  CV Reduction       : %.1f%% (median)\n", 
            median(power_improvement$CV_Reduction_Percent, na.rm = TRUE)))
cat(sprintf("  Power Gain (1.5x)  : +%.1f%% points\n", 
            mean(power_improvement$Power_Gain[power_improvement$Effect_Size == 1.5], na.rm = TRUE) * 100))
cat(sprintf("  Power Gain (2.0x)  : +%.1f%% points\n\n", 
            mean(power_improvement$Power_Gain[power_improvement$Effect_Size == 2.0], na.rm = TRUE) * 100))

# -----------------------------------------------------------------------------
# Step 6: Sample Size Recommendations
# -----------------------------------------------------------------------------

cat(rep("=", 80), "\n", sep = "")
cat("STEP 4: Sample Size Recommendations for Future Studies\n")
cat(rep("=", 80), "\n\n", sep = "")

sample_size_recs <- calculate_required_sample_size(
  power_stats = power_stats_corrected,
  effect_sizes = opt$effect_sizes,
  power_levels = opt$power_levels,
  alpha = opt$alpha
)

cat("Recommendations (based on corrected data):\n\n")
print(sample_size_recs[, c("Effect_Size", "Target_Power", "Required_N")], 
      row.names = FALSE)

cat(sprintf("\nBased on:\n"))
cat(sprintf("  Median CV    : %.3f\n", sample_size_recs$Median_CV[1]))
cat(sprintf("  Mean Depth   : %.1f\n\n", sample_size_recs$Mean_Depth[1]))

# -----------------------------------------------------------------------------
# Step 7: Save Results
# -----------------------------------------------------------------------------

cat(rep("=", 80), "\n", sep = "")
cat("💾 Saving results to disk...\n")
cat(rep("=", 80), "\n\n", sep = "")

# Save original power results
write.csv(
  actual_power_original,
  file.path(opt$outdir, "power_distribution_original.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  ✓ power_distribution_original.csv\n")

# Save corrected power results
write.csv(
  actual_power_corrected,
  file.path(opt$outdir, "power_distribution_corrected.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  ✓ power_distribution_corrected.csv\n")

# Save improvement comparison
write.csv(
  power_improvement,
  file.path(opt$outdir, "power_comparison_original_vs_corrected.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  ✓ power_comparison_original_vs_corrected.csv\n")

# Save combined dataset
power_combined <- bind_rows(actual_power_original, actual_power_corrected)
write.csv(
  power_combined,
  file.path(opt$outdir, "power_comparison_combined.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  ✓ power_comparison_combined.csv\n")

# Save sample size recommendations
write.csv(
  sample_size_recs,
  file.path(opt$outdir, "recommended_sample_sizes.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  ✓ recommended_sample_sizes.csv\n")

# Create summary statistics
summary_stats <- data.frame(
  Metric = c(
    "Total_Samples",
    "Analysis_Groups",
    "Biological_Conditions",
    "Genes_Analyzed_Original",
    "Genes_Analyzed_Corrected",
    "Median_CV_Original",
    "Median_CV_Corrected",
    "CV_Reduction_Percent",
    "Mean_Depth",
    "Mean_Power_Gain_1.5fold",
    "Mean_Power_Gain_2fold"
  ),
  Value = c(
    length(common_samples),
    length(groups),
    length(conditions),
    mean(power_stats_original$Genes_Analyzed, na.rm = TRUE),
    mean(power_stats_corrected$Genes_Analyzed, na.rm = TRUE),
    median(power_stats_original$CV, na.rm = TRUE),
    median(power_stats_corrected$CV, na.rm = TRUE),
    median(power_improvement$CV_Reduction_Percent, na.rm = TRUE),
    mean(power_stats_corrected$Mean_Depth, na.rm = TRUE),
    mean(power_improvement$Power_Gain[power_improvement$Effect_Size == 1.5], na.rm = TRUE),
    mean(power_improvement$Power_Gain[power_improvement$Effect_Size == 2.0], na.rm = TRUE)
  )
)

write.csv(
  summary_stats,
  file.path(opt$outdir, "power_comparison_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)
cat("  ✓ power_comparison_summary.csv\n")

# Create detailed log file
log_file <- file.path(opt$outdir, "power_analysis_log.txt")
sink(log_file)

cat("RNA-seq Power Analysis - Detailed Log\n")
cat(rep("=", 80), "\n", sep = "")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Experiment ID:", opt$experiment_id, "\n\n")

cat("CONFIGURATION:\n")
cat(rep("-", 80), "\n", sep = "")
cat("  Condition column     :", opt$condition_column, "\n")
cat("  Batch column         :", ifelse(is.null(opt$batch_column), "None", opt$batch_column), "\n")
cat("  Min count threshold  :", opt$min_count, "\n")
cat("  Min sample proportion:", opt$min_prop, "\n")
cat("  Effect sizes tested  :", paste(opt$effect_sizes, collapse = ", "), "\n")
cat("  Target power levels  :", paste(opt$power_levels, collapse = ", "), "\n")
cat("  Significance level   :", opt$alpha, "\n\n")

cat("INPUT DATA SUMMARY:\n")
cat(rep("-", 80), "\n", sep = "")
cat("  Total samples        :", length(common_samples), "\n")
cat("  Analysis groups      :", length(groups), "\n")
cat("  Biological conditions:", length(conditions), "\n")
cat("  Conditions           :", paste(conditions, collapse = ", "), "\n\n")

cat("ORIGINAL DATA (Before Batch Correction):\n")
cat(rep("-", 80), "\n", sep = "")
print(power_stats_original)
cat("\n")

cat("CORRECTED DATA (After Batch Correction):\n")
cat(rep("-", 80), "\n", sep = "")
print(power_stats_corrected)
cat("\n")

cat("POWER IMPROVEMENT:\n")
cat(rep("-", 80), "\n", sep = "")
print(power_improvement)
cat("\n")

cat("SAMPLE SIZE RECOMMENDATIONS:\n")
cat(rep("-", 80), "\n", sep = "")
print(sample_size_recs)
cat("\n")

cat("INTERPRETATION GUIDE:\n")
cat(rep("-", 80), "\n", sep = "")
cat("\nCV (Coefficient of Variation):\n")
cat("  - Measures within-group variability\n")
cat("  - Lower CV = more consistent measurements\n")
cat("  - Typical range: 0.1 (inbred animals) to 0.4 (human studies)\n\n")

cat("Power:\n")
cat("  - Probability of detecting true differences\n")
cat("  - 0.80 = 80% chance of detecting effect (minimum acceptable)\n")
cat("  - 0.90 = 90% chance of detecting effect (preferred)\n\n")

cat("Effect Size:\n")
cat("  - Fold change in gene expression\n")
cat("  - 1.5 = 50% increase (moderate biological effect)\n")
cat("  - 2.0 = 100% increase (large biological effect)\n\n")

cat("Sample Size Recommendations:\n")
cat("  - Required N is per group (multiply by number of groups for total)\n")
cat("  - Based on median CV and mean sequencing depth\n")
cat("  - Assumes similar experimental conditions\n\n")

cat(rep("=", 80), "\n", sep = "")
cat("Analysis completed successfully\n")
cat(rep("=", 80), "\n", sep = "")

sink()
cat("  ✓ power_analysis_log.txt\n\n")

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------

cat(rep("=", 80), "\n", sep = "")
cat("✅ Power Analysis Completed Successfully\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Key Findings:\n")
cat(sprintf("  • Batch correction reduced CV by %.1f%%\n", 
            median(power_improvement$CV_Reduction_Percent, na.rm = TRUE)))
cat(sprintf("  • Power improved by +%.1f%% points for 1.5-fold changes\n", 
            mean(power_improvement$Power_Gain[power_improvement$Effect_Size == 1.5], na.rm = TRUE) * 100))

# Handle case where no 1.5-fold data exists
if (any(sample_size_recs$Effect_Size == 1.5 & sample_size_recs$Target_Power == 0.8)) {
  recommended_n <- sample_size_recs$Required_N[
    sample_size_recs$Effect_Size == 1.5 & sample_size_recs$Target_Power == 0.8
  ][1]
  cat(sprintf("  • Recommended sample size: %d per group (for 1.5-fold, 80%% power)\n", recommended_n))
} else {
  cat("  • Sample size recommendations: See recommended_sample_sizes.csv\n")
}

cat("\nOutput Files:\n")
cat("  All results saved to:", opt$outdir, "\n")
cat("  Review 'power_analysis_log.txt' for detailed interpretation\n\n")

cat(rep("=", 80), "\n\n", sep = "")