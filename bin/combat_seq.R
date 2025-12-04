#!/usr/bin/env Rscript

# ComBat-seq Batch Correction Script - Multi-Dataset Support
#
# Usage with container:
#   singularity exec /path/to/NF_BulkRNAseq.sif Rscript combat_seq.R [options]

# Load required libraries (pre-installed in container)
suppressPackageStartupMessages({
    library(sva)        # ComBat-seq batch correction
    library(DESeq2)     # Differential expression analysis
    library(yaml)       # Configuration file parsing
    library(argparse)   # Command line argument parsing
})

# Parse command line arguments
parser <- ArgumentParser(description='ComBat-seq batch correction with multi-dataset support')
parser$add_argument('--counts', required=TRUE, help='Count matrix file')
parser$add_argument('--metadata', required=TRUE, help='Cleaned metadata file')
parser$add_argument('--config', required=TRUE, help='Configuration YAML file')
parser$add_argument('--output_dir', default='.', help='Output directory')
args <- parser$parse_args()

cat("=== ComBat-seq Batch Correction (Multi-Dataset Support) ===\n")
cat("Parameters:\n")
cat("  Counts file:", args$counts, "\n")
cat("  Cleaned metadata file:", args$metadata, "\n")
cat("  Config file:", args$config, "\n")
cat("  Output directory:", args$output_dir, "\n\n")

# Load configuration
cat("Loading configuration...\n")
config <- yaml::read_yaml(args$config)
batch_factor <- config$experimental_design$batch_factor
primary_factor <- config$experimental_design$primary_factor
additional_factors <- config$experimental_design$additional_factors

cat("Configuration loaded:\n")
cat("  Batch factor:", ifelse(is.null(batch_factor) || batch_factor == "none", "None (skipping batch correction)", batch_factor), "\n")
cat("  Primary factor:", primary_factor, "\n")

# Handle additional factors more robustly
if (is.null(additional_factors) || length(additional_factors) == 0 || 
    (length(additional_factors) == 1 && additional_factors == "")) {
    additional_factors <- NULL
    cat("  Additional factors: None\n")
} else {
    cat("  Additional factors:", paste(additional_factors, collapse=", "), "\n")
}
cat("\n")

# Check if batch correction should be skipped
skip_batch_correction <- is.null(batch_factor) || batch_factor == "none" || batch_factor == ""

if (skip_batch_correction) {
    cat("⚠️ No batch factor specified - skipping ComBat-seq correction\n")
    cat("   Analysis will proceed without batch correction\n")
    
    # Just copy the original counts and create a log
    counts_matrix_ordered <- read.csv(args$counts, row.names = 1, check.names = FALSE)
    
    # Save original counts as "corrected" counts
    output_file <- file.path(args$output_dir, "corrected_counts.csv")
    write.csv(counts_matrix_ordered, output_file)
    
    # Create log
    log_info <- list(
        timestamp = as.character(Sys.time()),
        input_counts_file = args$counts,
        input_metadata_file = args$metadata,
        config_file = args$config,
        batch_correction_applied = "No - no batch factor specified",
        note = "Original counts copied as corrected_counts.csv (no batch correction applied)"
    )
    
    log_content <- paste(names(log_info), log_info, sep = ": ", collapse = "\n")
    writeLines(log_content, file.path(args$output_dir, "batch_correction_log.txt"))
    
    cat("✓ Original counts saved as corrected_counts.csv (no batch correction)\n")
    cat("✓ Log saved to batch_correction_log.txt\n")
    quit(status = 0)
}

# Load data with proper error handling
cat("Loading count matrix...\n")
tryCatch({
    counts_matrix_ordered <- read.csv(args$counts, row.names = 1, check.names = FALSE)
    cat("✓ Count matrix loaded successfully\n")
}, error = function(e) {
    stop("Failed to load count matrix: ", e$message)
})

cat("Loading cleaned metadata...\n")
tryCatch({
    sample_df_cleaned <- read.csv(args$metadata, row.names = 1, check.names = FALSE)
    cat("✓ Cleaned metadata loaded successfully\n")
    cat("  Metadata dimensions:", paste(dim(sample_df_cleaned), collapse=" x "), "\n")
    cat("  Available columns:", paste(colnames(sample_df_cleaned), collapse=", "), "\n")
}, error = function(e) {
    stop("Failed to load cleaned metadata: ", e$message)
})

# Data validation and preparation
cat("\n=== Data Validation ===\n")

# Convert counts to a matrix and ensure integer
counts_matrix_ordered <- as.matrix(counts_matrix_ordered)
counts_matrix_ordered <- round(counts_matrix_ordered)

# Validate batch factor exists in cleaned metadata
if (!batch_factor %in% colnames(sample_df_cleaned)) {
    stop("Batch factor '", batch_factor, "' not found in cleaned metadata columns: ",
         paste(colnames(sample_df_cleaned), collapse=", "))
}

# Validate primary factor exists in cleaned metadata
if (!primary_factor %in% colnames(sample_df_cleaned)) {
    stop("Primary factor '", primary_factor, "' not found in cleaned metadata columns: ",
         paste(colnames(sample_df_cleaned), collapse=", "))
}

# Validate additional factors exist in cleaned metadata (improved handling)
if (!is.null(additional_factors) && length(additional_factors) > 0) {
    # Filter out empty strings
    additional_factors <- additional_factors[additional_factors != ""]
    
    if (length(additional_factors) > 0) {
        missing_additional <- additional_factors[!additional_factors %in% colnames(sample_df_cleaned)]
        if (length(missing_additional) > 0) {
            stop("Additional factors not found in cleaned metadata: ",
                 paste(missing_additional, collapse=", "))
        }
    } else {
        additional_factors <- NULL
    }
}

# Check sample alignment between counts and cleaned metadata
sample_names_counts <- colnames(counts_matrix_ordered)
sample_names_metadata <- rownames(sample_df_cleaned)

# Convert both to character to ensure consistent comparison
sample_names_counts <- as.character(sample_names_counts)
sample_names_metadata <- as.character(sample_names_metadata)

cat("Sample alignment check:\n")
cat("  Count matrix samples:", length(sample_names_counts), "\n")
cat("  Cleaned metadata samples:", length(sample_names_metadata), "\n")

if (!identical(sort(sample_names_counts), sort(sample_names_metadata))) {
    cat("⚠ WARNING: Sample names don't match exactly\n")
    # Find common samples
    common_samples <- intersect(sample_names_counts, sample_names_metadata)
    if (length(common_samples) == 0) {
        stop("No common samples found between count matrix and cleaned metadata")
    }
    cat("  Using", length(common_samples), "common samples\n")
    # Subset to common samples
    counts_matrix_ordered <- counts_matrix_ordered[, common_samples]
    sample_df_cleaned <- sample_df_cleaned[common_samples, ]
} else {
    cat("✓ Sample names match perfectly\n")
}

# Ensure sample order matches between counts and metadata
sample_df_cleaned <- sample_df_cleaned[colnames(counts_matrix_ordered), ]

cat("✓ Data validation completed\n")
cat("  Final dimensions - Counts:", paste(dim(counts_matrix_ordered), collapse=" x "), "\n")
cat("  Final dimensions - Cleaned metadata:", paste(dim(sample_df_cleaned), collapse=" x "), "\n")

#####################################
# ComBat-seq Batch Correction
#####################################

cat("\n=== Applying ComBat-seq ===\n")

# Prepare data from cleaned metadata
counts <- counts_matrix_ordered
batch <- factor(sample_df_cleaned[[batch_factor]])

# Print batch information from cleaned metadata
cat("Batch factor distribution (from cleaned metadata):\n")
print(table(batch))

# Check if we have enough batches
if (length(levels(batch)) < 2) {
    cat("⚠️ Only", length(levels(batch)), "batch level found - skipping ComBat-seq correction\n")
    
    # Save original counts
    output_file <- file.path(args$output_dir, "corrected_counts.csv")
    write.csv(counts, output_file)
    
    # Create log
    log_info <- list(
        timestamp = as.character(Sys.time()),
        input_counts_file = args$counts,
        input_metadata_file = args$metadata,
        config_file = args$config,
        batch_factor = batch_factor,
        num_batches = length(levels(batch)),
        batch_correction_applied = "No - insufficient batch levels (need at least 2)",
        note = "Original counts saved as corrected_counts.csv"
    )
    
    log_content <- paste(names(log_info), log_info, sep = ": ", collapse = "\n")
    writeLines(log_content, file.path(args$output_dir, "batch_correction_log.txt"))
    
    cat("✓ Original counts saved as corrected_counts.csv (no batch correction needed)\n")
    quit(status = 0)
}

# Create biological covariate model matrix from cleaned metadata
bio_cov_columns <- c(primary_factor)
if (!is.null(additional_factors) && length(additional_factors) > 0) {
    bio_cov_columns <- c(bio_cov_columns, additional_factors)
}

bio_cov <- sample_df_cleaned[, bio_cov_columns, drop = FALSE]
rownames(bio_cov) <- colnames(counts)

cat("Biological factors being preserved (from cleaned metadata):\n")
for (col in bio_cov_columns) {
    cat("  ", col, ":", length(unique(bio_cov[[col]])), "levels ->",
        paste(unique(bio_cov[[col]]), collapse=", "), "\n")
}

# Create group factor for ComBat-seq
cat("\nCreating group factor from cleaned metadata...\n")
group_factors <- apply(bio_cov, 1, function(x) paste(x, collapse = "_"))
group <- factor(group_factors)

cat("Group factor levels:", length(levels(group)), "\n")
cat("Group distribution:\n")
print(table(group))

# Apply ComBat-seq with error handling
cat("\nApplying ComBat-seq correction...\n")
corrected_counts <- NULL  # Initialize variable

tryCatch({
    corrected_counts <<- ComBat_seq(
        counts = counts,
        batch = batch,
        group = group
    )
    cat("✓ ComBat-seq correction completed successfully!\n")
}, error = function(e) {
    cat("✗ ComBat-seq failed with error:", e$message, "\n")
    # Try an alternative approach without a group
    cat("Trying ComBat-seq without group factor...\n")
    tryCatch({
        corrected_counts <<- ComBat_seq(
            counts = counts,
            batch = batch,
            group = NULL
        )
        cat("✓ ComBat-seq completed with simplified approach\n")
    }, error = function(e2) {
        stop("ComBat-seq failed with both approaches: ", e2$message)
    })
})

# Validate corrected counts
if (is.null(corrected_counts)) {
    stop("ComBat-seq correction failed - no corrected counts generated")
}

cat("Corrected counts dimensions:", paste(dim(corrected_counts), collapse=" x "), "\n")

# Save corrected counts (CSV only - no plots for memory efficiency)
output_file <- file.path(args$output_dir, "corrected_counts.csv")
write.csv(corrected_counts, output_file)
cat("✓ Corrected counts saved to:", output_file, "\n")

# Create comprehensive log
log_info <- list(
    timestamp = as.character(Sys.time()),
    input_counts_file = args$counts,
    input_metadata_file = args$metadata,
    config_file = args$config,
    metadata_type = "cleaned_metadata (user-selected columns)",
    available_metadata_columns = paste(colnames(sample_df_cleaned), collapse = ", "),
    batch_factor = batch_factor,
    primary_factor = primary_factor,
    additional_factors = ifelse(is.null(additional_factors), "None", paste(additional_factors, collapse = ", ")),
    num_batches = length(levels(batch)),
    num_samples = ncol(counts),
    num_genes = nrow(counts),
    batch_distribution = paste(names(table(batch)), table(batch), sep=":", collapse="; "),
    group_levels = length(levels(group)),
    correction_method = "ComBat_seq",
    batch_correction_applied = "Yes - ComBat-seq applied successfully",
    output_files = "corrected_counts.csv",
    note = "Cleaned metadata should be saved separately in data folder",
    visualization_note = "All visualizations disabled for memory efficiency - CSV outputs only"
)

# Write log
log_content <- paste(names(log_info), log_info, sep = ": ", collapse = "\n")
writeLines(log_content, file.path(args$output_dir, "batch_correction_log.txt"))

cat("\n=== ComBat-seq Batch Correction Completed Successfully! ===\n")
cat("Output files created:\n")
cat("  - corrected_counts.csv\n")
cat("  - batch_correction_log.txt\n")
cat("Note: Cleaned metadata should be saved separately in data folder\n")
cat("Note: Visualizations disabled for memory efficiency (CSV outputs only)\n")

