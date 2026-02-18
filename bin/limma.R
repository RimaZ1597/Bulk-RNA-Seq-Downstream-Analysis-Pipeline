#!/usr/bin/env Rscript

# limma.R
# Differential Expression Analysis using limma-voom
# Outputs: ALL genes (no p-value filtering) + normalized counts + contrast-specific results

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(yaml)
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# Parse command line arguments
parser <- ArgumentParser(description='limma-voom differential expression analysis')
parser$add_argument('--counts', required=TRUE,
                    help='Path to corrected counts CSV (from ComBat-seq or filtered counts)')
parser$add_argument('--metadata', required=TRUE,
                    help='Path to sample metadata CSV')
parser$add_argument('--config', required=TRUE,
                    help='Path to analysis_config.yaml (from setup.py)')
parser$add_argument('--output_dir', default='.',
                    help='Output directory for results (NextFlow publishDir handles path structure)')
parser$add_argument('--contrast', default='all',
                    help='Specific contrast name to run (or "all" for all contrasts)')
parser$add_argument('--combatseq_applied', default='false',
                    help='Flag indicating if ComBat-seq was applied (true/false)')
args <- parser$parse_args()

cat("\n=== limma-voom Differential Expression Analysis ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Input Parameters:\n")
cat("  Counts matrix:", args$counts, "\n")
cat("  Sample metadata:", args$metadata, "\n")
cat("  Configuration:", args$config, "\n")
cat("  Output directory:", args$output_dir, "\n")
cat("  Contrast filter:", args$contrast, "\n")
cat("  ComBat-seq applied:", args$combatseq_applied, "\n\n")

# Load configuration
cat("Loading analysis configuration...\n")
tryCatch({
  config <- yaml::read_yaml(args$config)
}, error = function(e) {
  stop("Failed to load configuration file: ", e$message)
})

# Extract parameters
primary_factor <- config$experimental_design$primary_factor
batch_factor <- config$experimental_design$batch_factor
additional_factors <- config$experimental_design$additional_factors
project_name <- config$project_name
comparisons <- config$comparisons

cat("Experimental Design Configuration:\n")
cat("  Project:", project_name, "\n")
cat("  Primary factor:", primary_factor, "\n")
cat("  Batch factor:", ifelse(is.null(batch_factor), "None", batch_factor), "\n")
cat("  Additional factors:", ifelse(length(additional_factors) > 0,
                                    paste(additional_factors, collapse=", "), "None"), "\n")
cat("  Total contrasts defined:", length(comparisons), "\n\n")

# Validate DE is requested
if (!config$analysis_parameters$perform_differential_expression) {
  stop("Differential expression not requested in configuration.")
}

# Natural sorting function for alphanumeric column names
natural_sort <- function(x) {
  # Extract numeric parts and sort naturally
  x[order(nchar(x), x)]
}

# Advanced natural sort that handles mixed alphanumeric correctly
mixedsort_base <- function(x) {
  # Split each string into parts (letters and numbers)
  split_parts <- strsplit(x, "(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)", perl = TRUE)
  
  # Create sorting keys
  keys <- lapply(split_parts, function(parts) {
    # Convert numeric parts to integers, keep text parts as-is
    sapply(parts, function(part) {
      if (grepl("^\\d+$", part)) {
        as.numeric(part)
      } else {
        part
      }
    }, USE.NAMES = FALSE)
  })
  
  # Sort using the keys
  order_idx <- order(sapply(keys, function(k) paste(sprintf("%010s", k), collapse = "")))
  x[order_idx]
}

# Load count matrix (CSV: genes as rows, samples as columns)
cat("Loading count matrix...\n")
tryCatch({
  counts_df <- read.csv(args$counts, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  counts_df <- as.data.frame(lapply(counts_df, function(x) as.numeric(as.character(x))))
  rownames(counts_df) <- rownames(read.csv(args$counts, row.names = 1, check.names = FALSE))
}, error = function(e) {
  stop("Failed to load counts matrix: ", e$message)
})

cat("  Dimensions (before sorting):", nrow(counts_df), "genes x", ncol(counts_df), "samples\n")

# Apply natural sorting to column names to fix P1wellA1, P1wellA10, P1wellA2 -> P1wellA1, P1wellA2, P1wellA10
cat("Applying natural sort to column names...\n")
original_order <- colnames(counts_df)
natural_order <- mixedsort_base(original_order)
counts_df <- counts_df[, natural_order]

cat("  Column order fixed: ", paste(head(colnames(counts_df), 3), collapse = ", "), "...\n")
cat("  Dimensions (after sorting):", nrow(counts_df), "genes x", ncol(counts_df), "samples\n")

# Load metadata (CSV)
cat("Loading sample metadata...\n")
tryCatch({
  metadata_df <- read.csv(args$metadata, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Apply the same natural sorting to metadata to keep it aligned with counts
  cat("Aligning metadata with naturally sorted column order...\n")
  metadata_df <- metadata_df[natural_order, , drop = FALSE]
  
}, error = function(e) {
  stop("Failed to load metadata: ", e$message)
})

cat("  Samples in metadata:", nrow(metadata_df), "\n")

# Verify that counts columns and metadata rows are perfectly aligned
if (!all(colnames(counts_df) == rownames(metadata_df))) {
  stop("CRITICAL ERROR: Column names in counts matrix do not match row names in metadata after sorting!")
} else {
  cat("✓ Counts matrix columns and metadata rows are perfectly aligned\n")
}

cat("\n")

# === CRITICAL: Sample Order Validation ===
cat("Validating sample consistency...\n")
counts_samples <- colnames(counts_df)
metadata_samples <- rownames(metadata_df)

if (!all(counts_samples %in% metadata_samples)) {
  missing <- setdiff(counts_samples, metadata_samples)
  stop("Samples in counts not found in metadata: ", paste(missing, collapse=", "))
}

# Align samples (same order as counts matrix)
common_samples <- intersect(counts_samples, metadata_samples)
counts_matrix <- as.matrix(counts_df[, common_samples])
metadata_aligned <- metadata_df[common_samples, ]

cat("  Final aligned dataset:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "samples\n")
cat("  Sample order validated:", identical(colnames(counts_matrix), rownames(metadata_aligned)), "\n\n")

if (!identical(colnames(counts_matrix), rownames(metadata_aligned))) {
  stop("CRITICAL: Sample order mismatch!")
}

# === QUALITY CHECKS ===
cat("=== Quality Checks ===\n")
cat("Zeros:", sum(counts_matrix == 0),
    "(", round(100*sum(counts_matrix == 0)/length(counts_matrix), 2), "%)\n")
cat("Negative values:", sum(counts_matrix < 0), "\n")
cat("NA values:", sum(is.na(counts_matrix)), "\n")
cat("Library sizes - Min:", min(colSums(counts_matrix)),
    "Max:", max(colSums(counts_matrix)), "\n\n")

# Validate required factors
required_factors <- c(primary_factor)
if (!is.null(batch_factor) && batch_factor != "None") {
  required_factors <- c(required_factors, batch_factor)
}

missing_factors <- setdiff(required_factors, colnames(metadata_aligned))
if (length(missing_factors) > 0) {
  stop("Required factors not in metadata: ", paste(missing_factors, collapse=", "))
}

# === CREATE DGEList AND NORMALIZE ===
cat("=== Preparing Data for limma-voom ===\n")

dge <- DGEList(counts = counts_matrix,
               samples = metadata_aligned,
               genes = data.frame(gene_symbol = rownames(counts_matrix), stringsAsFactors = FALSE))

cat("DGEList created: ", nrow(dge), "genes x", ncol(dge), "samples\n")

# Calculate TMM normalization factors
cat("Calculating TMM normalization factors...\n")
dge <- calcNormFactors(dge, method = "TMM")

cat("TMM normalization factors:\n")
cat("  Range:", round(min(dge$samples$norm.factors), 3), "to",
    round(max(dge$samples$norm.factors), 3), "\n")
cat("  Mean:", round(mean(dge$samples$norm.factors), 3), "\n\n")

# === SAVE NORMALIZED COUNTS ===
cat("=== Saving Normalized Count Matrices ===\n")

# 1. TMM-normalized CPM (non-logged, for visualization/input to other tools)
cpm_normalized <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)
cpm_df <- data.frame(
  gene_symbol = rownames(cpm_normalized),
  cpm_normalized,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(cpm_df,
          file.path(args$output_dir, "normalized_counts_cpm.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Saved: normalized_counts_cpm.csv (TMM-normalized CPM, non-logged)\n")

# 2. Log2-CPM (voom will produce this, save for reference)
log2cpm_normalized <- cpm(dge, log = TRUE, prior.count = 2)
log2cpm_df <- data.frame(
  gene_symbol = rownames(log2cpm_normalized),
  log2cpm_normalized,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(log2cpm_df,
          file.path(args$output_dir, "normalized_counts_log2cpm.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Saved: normalized_counts_log2cpm.csv (log2-CPM for reference)\n\n")

# === BUILD UNIVERSAL MULTI-FACTOR DESIGN MATRIX ===
cat("=== Setting Up Multi-factor Design ===\n")

# UNIVERSAL APPROACH: Automatically detect all factors used in contrasts
contrast_factors <- c()
for (contrast in comparisons) {
  # Extract all factor names from contrast format like "factor1:value1_factor2:value2_vs_value3"
  # This regex finds all "word:" patterns
  factor_matches <- regmatches(contrast, gregexpr("\\b[a-zA-Z_][a-zA-Z0-9_]*(?=:)", contrast, perl = TRUE))[[1]]
  if (length(factor_matches) > 0) {
    contrast_factors <- c(contrast_factors, factor_matches)
  }
}
contrast_factors <- unique(contrast_factors)

cat("Factors detected in contrasts:", paste(contrast_factors, collapse=", "), "\n")

# Build design factor list: primary factor + detected contrast factors
design_factors <- unique(c(primary_factor, contrast_factors))

# Validate all factors exist in metadata
available_factors <- intersect(design_factors, colnames(metadata_aligned))
missing_factors <- setdiff(design_factors, colnames(metadata_aligned))

if (length(missing_factors) > 0) {
  cat("⚠ Warning: Factors in contrasts not found in metadata:", paste(missing_factors, collapse=", "), "\n")
  cat("  Available factors:", paste(colnames(metadata_aligned), collapse=", "), "\n")
  design_factors <- available_factors
}

cat("Final design factors:", paste(design_factors, collapse=" + "), "\n")

# Create multi-factor treatment combinations
if (length(design_factors) == 1) {
  # Single factor design
  Treat <- factor(metadata_aligned[[design_factors[1]]])
  cat("✓ Single-factor design:", design_factors[1], "\n")
} else {
  # Multi-factor design: combine all factors with "." separator
  combined_values <- apply(metadata_aligned[, design_factors, drop=FALSE], 1, 
                          function(x) paste(x, collapse="."))
  Treat <- factor(combined_values)
  cat("✓ Multi-factor design:", paste(design_factors, collapse=" × "), "\n")
}

cat("Treatment factor levels:\n")
treat_table <- table(Treat)
print(treat_table)

# Create design matrix (no intercept)
design <- model.matrix(~0 + Treat)

# Clean column names
original_names <- levels(Treat)
clean_names <- make.names(levels(Treat))
colnames(design) <- clean_names

cat("\nDesign matrix: ", nrow(design), "x", ncol(design), "\n")
cat("Columns (first 5):", paste(head(clean_names, 5), collapse=", "), "\n\n")

# Name mapping for contrast parsing
name_mapping <- data.frame(
  Original = original_names,
  Clean = clean_names,
  stringsAsFactors = FALSE
)

# === DESIGN MATRIX DIAGNOSTIC (MOVED HERE - BEFORE PARSING) ===
print_design_matrix_info <- function(name_mapping) {
  cat("\n=== Design Matrix Column Information ===\n")
  cat("Total columns:", nrow(name_mapping), "\n\n")
  cat("Sample of columns (showing original -> cleaned names):\n")
  for (i in seq_len(min(10, nrow(name_mapping)))) {
    cat(sprintf("  %2d. Original: %-50s\n", i, name_mapping$Original[i]))
    cat(sprintf("      Clean:    %-50s\n", name_mapping$Clean[i]))
    cat("\n")
  }
  if (nrow(name_mapping) > 10) {
    cat("  ... (", nrow(name_mapping) - 10, "more columns)\n\n")
  }
}

print_design_matrix_info(name_mapping)

# === VOOM TRANSFORMATION ===
cat("=== Applying voom Transformation ===\n")
v <- voom(dge, design = design, plot = FALSE)
cat("✓ Voom transformation completed\n")
cat("  Expression matrix:", dim(v$E), "\n\n")

# Save voom plot
png(file.path(args$output_dir, "voom_mean_variance_plot.png"),
    width = 10, height = 6, units = "in", res = 300)
v_plot <- voom(dge, design = design, plot = TRUE)
dev.off()
cat("✓ Voom plot saved\n")

# Save voom log2-CPM values (this is what limma actually uses)
voom_log2cpm_df <- data.frame(
  gene_symbol = rownames(v$E),
  v$E,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.csv(voom_log2cpm_df,
          file.path(args$output_dir, "voom_log2cpm.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Saved: voom_log2cpm.csv (voom-transformed expression used in limma)\n\n")

# === ESTIMATE WITHIN-SUBJECT CORRELATION ===
use_blocking <- FALSE
corfit <- NULL

if (!is.null(batch_factor) && batch_factor != "None" && batch_factor %in% colnames(metadata_aligned)) {
  cat("=== Estimating Within-Subject Correlation ===\n")
  cat("Blocking factor:", batch_factor, "\n")
  block_vector <- factor(metadata_aligned[[batch_factor]])
  tryCatch({
    corfit <- duplicateCorrelation(v, design, block = block_vector)
    cat("Consensus correlation:", round(corfit$consensus, 4), "\n")
    if (corfit$consensus > 0.01) {
      use_blocking <- TRUE
      cat("✓ Will use blocking in linear model\n\n")
    } else {
      cat("ℹ Correlation too low - no blocking\n\n")
    }
  }, error = function(e) {
    cat("⚠ Could not estimate correlation:", e$message, "\n")
    cat("  Proceeding without blocking\n\n")
  })
}

# === FIT LINEAR MODEL ===
cat("=== Fitting Linear Model ===\n")

if (use_blocking && !is.null(corfit)) {
  fit <- lmFit(v, design,
               block = factor(metadata_aligned[[batch_factor]]),
               correlation = corfit$consensus)
  cat("✓ Model fitted WITH blocking\n")
} else {
  fit <- lmFit(v, design)
  cat("✓ Model fitted WITHOUT blocking\n")
}

cat("  Genes:", nrow(fit$coefficients), "\n")
cat("  Coefficients:", ncol(fit$coefficients), "\n\n")

# === UNIVERSAL CONTRAST PARSING FOR ANY DATASET ===
parse_contrast_robust <- function(contrast_str, name_mapping) {
  # Universal parser for any contrast format:
  # "factor1:value1_factor2:value2_vs_factor3:value3"
  # "factor1:value1_vs_value2"  (within-factor comparison)
  # "factor1:value1_factor2:value2_vs_value3_value4"  (mixed format)
  
  cat("\n--- Parsing contrast:", contrast_str, "\n")
  
  # Split by "_vs_" to separate the two groups
  parts <- strsplit(contrast_str, "_vs_")[[1]]
  if (length(parts) != 2) {
    warning("Cannot parse contrast (missing _vs_): ", contrast_str)
    return(NULL)
  }
  
  group1_str <- parts[1]
  group2_str <- parts[2]
  
  cat("  Group 1 string:", group1_str, "\n")
  cat("  Group 2 string:", group2_str, "\n")
  
  # === UNIVERSAL FACTOR:VALUE EXTRACTION ===
  extract_factor_values <- function(group_str) {
    factor_values <- list()
    
    # Extract all factor:value pairs
    # Pattern: word followed by colon, then capture until next underscore or end
    pattern <- "([a-zA-Z_][a-zA-Z0-9_]*):([^_]+?)(?=_[a-zA-Z_][a-zA-Z0-9_]*:|$)"
    matches <- gregexpr(pattern, group_str, perl = TRUE)
    
    if (matches[[1]][1] != -1) {
      match_data <- regmatches(group_str, matches)[[1]]
      
      for (match in match_data) {
        # Split on first colon only
        colon_pos <- regexpr(":", match)
        if (colon_pos > 0) {
          factor_name <- substr(match, 1, colon_pos - 1)
          factor_value <- substr(match, colon_pos + 1, nchar(match))
          # Clean trailing underscore if present
          factor_value <- sub("_$", "", factor_value)
          factor_values[[factor_name]] <- factor_value
        }
      }
    }
    
    # Handle case where no factor:value pairs found (just a value)
    if (length(factor_values) == 0) {
      # Check if this looks like a standalone value (no colons)
      if (!grepl(":", group_str)) {
        cat("    ℹ Treating as standalone value:", group_str, "\n")
        return(list(standalone_value = group_str))
      }
    }
    
    return(factor_values)
  }
  
  # Extract factor values for both groups
  group1_factors <- extract_factor_values(group1_str)
  group2_factors <- extract_factor_values(group2_str)
  
  cat("  Group 1 factors:", 
      if(length(group1_factors) > 0) paste(names(group1_factors), "=", unlist(group1_factors), collapse=", ") else "none", "\n")
  cat("  Group 2 factors:", 
      if(length(group2_factors) > 0) paste(names(group2_factors), "=", unlist(group2_factors), collapse=", ") else "none", "\n")
  
  # === HANDLE WITHIN-FACTOR CONTRASTS ===
  # If group2 has no factor:value pairs, inherit factors from group1
  if (length(group2_factors) == 1 && "standalone_value" %in% names(group2_factors)) {
    if (length(group1_factors) > 0) {
      cat("    ℹ Within-factor contrast detected - inheriting factors from group 1\n")
      
      # Store the standalone value before copying
      standalone_val <- group2_factors[["standalone_value"]]
      
      # Copy all factors from group1, but replace the last factor's value
      group2_factors <- group1_factors
      factor_names <- names(group1_factors)
      
      # Replace the last factor's value with the standalone value
      last_factor <- factor_names[length(factor_names)]
      group2_factors[[last_factor]] <- standalone_val
      
      cat("    ℹ Updated group 2 factors:", 
          paste(names(group2_factors), "=", unlist(group2_factors), collapse=", "), "\n")
    }
  }
  
  # === BUILD EXPECTED DESIGN MATRIX COLUMN NAMES ===
  build_column_name <- function(factor_list) {
    if (length(factor_list) == 0) return(NULL)
    
    # Sort factors by name for consistency
    sorted_factors <- factor_list[order(names(factor_list))]
    values <- unlist(sorted_factors)
    
    return(paste(values, collapse = "."))
  }
  
  expected_col1 <- build_column_name(group1_factors)
  expected_col2 <- build_column_name(group2_factors)
  
  if (is.null(expected_col1) || is.null(expected_col2)) {
    cat("    ✗ Could not build expected column names\n")
    return(NULL)
  }
  
  cat("  Looking for design columns:", expected_col1, "vs", expected_col2, "\n")
  
  # === UNIVERSAL COLUMN MATCHING ===
  find_column_universal <- function(expected_name, orig_names, clean_names) {
    # Try exact match first
    exact_idx <- which(orig_names == expected_name)
    if (length(exact_idx) > 0) {
      cat("    ✓ Found exact match:", clean_names[exact_idx[1]], "\n")
      return(clean_names[exact_idx[1]])
    }
    
    # Try fuzzy match with regex escaping
    escaped_name <- gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", expected_name)
    fuzzy_matches <- grep(escaped_name, orig_names, ignore.case = TRUE, fixed = FALSE)
    
    if (length(fuzzy_matches) > 0) {
      cat("    ⚠ Found fuzzy match:", clean_names[fuzzy_matches[1]], "\n")
      return(clean_names[fuzzy_matches[1]])
    }
    
    # Try partial matching on individual components
    name_parts <- strsplit(expected_name, "\\.")[[1]]
    for (i in seq_along(orig_names)) {
      orig_name <- orig_names[i]
      if (all(sapply(name_parts, function(part) grepl(part, orig_name, fixed = TRUE)))) {
        cat("    ⚠ Found component match:", clean_names[i], "\n")
        return(clean_names[i])
      }
    }
    
    cat("    ✗ No match found for:", expected_name, "\n")
    return(NULL)
  }
  
  # Find the two columns
  col1_clean <- find_column_universal(expected_col1, name_mapping$Original, name_mapping$Clean)
  col2_clean <- find_column_universal(expected_col2, name_mapping$Original, name_mapping$Clean)
  
  if (is.null(col1_clean) || is.null(col2_clean)) {
    warning("Failed to map contrast to design matrix: ", contrast_str)
    cat("\n  Available design columns (first 10):\n")
    print(head(name_mapping, 10))
    return(NULL)
  }
  
  # === BUILD CONTRAST EXPRESSION ===
  contrast_expr <- paste(col1_clean, "-", col2_clean)
  cat("  ✓ Mapped to contrast expression:", contrast_expr, "\n")
  
  return(list(
    expression = contrast_expr,
    group1_cols = col1_clean,
    group2_cols = col2_clean,
    original_contrast = contrast_str
  ))
}

# === PARSE AND BUILD CONTRASTS ===
cat("=== Processing Contrasts ===\n")

# Filter contrasts if specific one requested
if (args$contrast != "all") {
  comparisons <- comparisons[comparisons == args$contrast]
  if (length(comparisons) == 0) {
    stop("Requested contrast not found: ", args$contrast)
  }
  cat("Running single contrast:", args$contrast, "\n")
}

# Parse all contrasts with detailed logging
parsed_contrasts <- list()
failed_contrasts <- c()

for (contrast_name in comparisons) {
  parsed <- parse_contrast_robust(contrast_name, name_mapping)
  if (!is.null(parsed)) {
    parsed_contrasts[[contrast_name]] <- parsed
    cat("✓ Successfully parsed:", contrast_name, "\n")
  } else {
    failed_contrasts <- c(failed_contrasts, contrast_name)
    cat("✗ Failed to parse:", contrast_name, "\n")
  }
}

# Report parsing results
cat("\n=== Contrast Parsing Summary ===\n")
cat("Total contrasts requested:", length(comparisons), "\n")
cat("Successfully parsed:", length(parsed_contrasts), "\n")
cat("Failed to parse:", length(failed_contrasts), "\n")

if (length(failed_contrasts) > 0) {
  cat("\nFailed contrasts:\n")
  for (fc in failed_contrasts) {
    cat("  -", fc, "\n")
  }
  cat("\nAvailable design matrix columns (first 10):\n")
  print(head(name_mapping, 10))
}

if (length(parsed_contrasts) == 0) {
  stop("No valid contrasts could be parsed! Check contrast names and design matrix.")
}

cat("\nProceeding with", length(parsed_contrasts), "valid contrasts\n\n")

# Build contrast matrix using only successfully parsed contrasts
contrast_expressions <- sapply(parsed_contrasts, function(x) x$expression)
contrast_names <- names(parsed_contrasts)

cat("=== Building Contrast Matrix ===\n")
cat("Contrast expressions:\n")
for (i in seq_along(contrast_names)) {
  cat(sprintf("  %d. %s = %s\n", i, contrast_names[i], contrast_expressions[i]))
}

# Create makeContrasts call
contrast_call <- paste0(
  "makeContrasts(",
  paste(paste0("`", contrast_names, "` = ", contrast_expressions), collapse = ", "),
  ", levels = design)"
)

cat("\nExecuting makeContrasts...\n")
tryCatch({
  cm <- eval(parse(text = contrast_call))
  cat("✓ Contrast matrix created successfully\n")
  cat("  Dimensions:", nrow(cm), "x", ncol(cm), "\n")
  cat("  Contrasts:\n")
  if (ncol(cm) > 0) {
    for (i in seq_len(ncol(cm))) {
      cat(sprintf("    %d. %s\n", i, colnames(cm)[i]))
    }
  }
  cat("\n")
}, error = function(e) {
  cat("ERROR creating contrast matrix:\n")
  cat("  Message:", e$message, "\n")
  cat("  Call attempted:", contrast_call, "\n")
  cat("\nAvailable design columns:\n")
  print(colnames(design))
  stop("Failed to create contrast matrix. Check contrast definitions.")
})

# === FIT CONTRASTS ===
cat("=== Fitting Contrasts ===\n")

fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

cat("✓ Empirical Bayes statistics computed\n")
cat("  Genes:", nrow(fit2), "\n")
cat("  Contrasts:", ncol(fit2), "\n\n")

# === SAVE RESULTS ===
cat("=== Saving Results ===\n")

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Thresholds (for filtering significant genes subset only)
fdr_threshold <- 0.05
logfc_threshold <- 1

# === SAVE GENERAL RESULTS (LIMMA MAIN FOLDER) ===
cat("\n=== Saving General Results to Main Limma Folder ===\n")

# Save normalized counts (log2-CPM from voom)
norm_counts <- v$E
rownames(norm_counts) <- rownames(counts_df)
colnames(norm_counts) <- colnames(counts_df)
write.csv(norm_counts,
          file.path(args$output_dir, "normalized_counts_log2cpm.csv"),
          quote = FALSE)
cat("✓ Normalized counts saved: normalized_counts_log2cpm.csv\n")

# Save design matrix information
write.csv(name_mapping,
          file.path(args$output_dir, "design_matrix_mapping.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Design matrix mapping saved\n")

# Initialize summary
de_summary <- data.frame(
  Contrast = character(),
  Total_Genes = integer(),
  Significant_FDR05_FC1 = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  stringsAsFactors = FALSE
)

# Store all results for combined output
all_results_combined <- list()

# === PROCESS CONTRAST-SPECIFIC RESULTS ===
cat("\n=== Processing Contrast-Specific Results ===\n")

# Process each contrast
if (ncol(cm) > 0) {
  for (i in seq_len(ncol(cm))) {
    contrast_name <- colnames(cm)[i]
    # Create clean folder name from contrast
    clean_folder_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)
    clean_folder_name <- gsub("_{2,}", "_", clean_folder_name)  # Remove multiple underscores
    clean_folder_name <- gsub("^_|_$", "", clean_folder_name)   # Remove leading/trailing underscores
    
    cat("\n--- Processing Contrast:", contrast_name, "---\n")
    cat("    Folder name:", clean_folder_name, "\n")
    
    # Create contrast-specific subdirectory
    contrast_dir <- file.path(args$output_dir, clean_folder_name)
    dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
    cat("    Results folder:", contrast_dir, "\n")
    
    # === GET ALL GENES (NO P-VALUE FILTERING) ===
    all_results <- topTable(fit2, coef = i, number = Inf, sort.by = "none")
    all_results$gene <- rownames(all_results)
    all_results <- all_results[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    
    # Add contrast column for combined output
    all_results$contrast <- contrast_name
    all_results_combined[[contrast_name]] <- all_results
    
    # === SAVE CONTRAST-SPECIFIC RESULTS ===
    # 1. All genes (complete results)
    write.csv(all_results,
              file.path(contrast_dir, "all_genes_results.csv"),
              row.names = FALSE, quote = FALSE)
    cat("    ✓ All genes saved:", nrow(all_results), "genes\n")
    
    # 2. Top genes by statistical significance
    n_top <- min(100, nrow(all_results))
    if (n_top > 0) {
      top100 <- all_results[order(all_results$P.Value), ][seq_len(n_top), ]
      write.csv(top100,
                file.path(contrast_dir, "top100_by_pvalue.csv"),
                row.names = FALSE, quote = FALSE)
      cat("    ✓ Top 100 genes by p-value saved\n")
    }
    
    # === SAVE SIGNIFICANT SUBSET (for convenience) ===
    sig_results <- all_results[all_results$adj.P.Val < fdr_threshold &
                                abs(all_results$logFC) > logfc_threshold, ]
    n_up <- sum(all_results$adj.P.Val < fdr_threshold & all_results$logFC > logfc_threshold)
    n_down <- sum(all_results$adj.P.Val < fdr_threshold & all_results$logFC < -logfc_threshold)
    n_sig <- n_up + n_down
    
    # 3. Significant genes (filtered results)
    if (nrow(sig_results) > 0) {
      write.csv(sig_results,
                file.path(contrast_dir, "significant_genes_FDR05_FC1.csv"),
                row.names = FALSE, quote = FALSE)
      cat("    ✓ Significant genes:", n_sig, "(", n_up, "up,", n_down, "down)\n")
    } else {
      cat("    ℹ No significant genes at FDR<0.05 & |logFC|>1\n")
    }
    
    # === CREATE VOLCANO PLOT ===
    all_results$significance <- "Not Significant"
    all_results$significance[all_results$adj.P.Val < fdr_threshold &
                              all_results$logFC > logfc_threshold] <- "Upregulated"
    all_results$significance[all_results$adj.P.Val < fdr_threshold &
                              all_results$logFC < -logfc_threshold] <- "Downregulated"
    
    all_results$neg_log10_fdr <- -log10(all_results$adj.P.Val)
    max_y <- max(all_results$neg_log10_fdr[!is.infinite(all_results$neg_log10_fdr)])
    all_results$neg_log10_fdr[is.infinite(all_results$neg_log10_fdr)] <- max_y + 10
    
    top_genes <- head(all_results[all_results$significance != "Not Significant", ], 10)
    
    p <- ggplot(all_results, aes(x = logFC, y = neg_log10_fdr, color = significance)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(
        values = c("Upregulated" = "#d62728",
                  "Downregulated" = "#1f77b4",
                  "Not Significant" = "grey70"),
        labels = c("Upregulated" = paste0("Upregulated (", n_up, ")"),
                  "Downregulated" = paste0("Downregulated (", n_down, ")"),
                  "Not Significant" = paste0("Not Significant (", nrow(all_results) - n_sig, ")"))
      ) +
      geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(fdr_threshold),
                linetype = "dashed", color = "grey40") +
      labs(
        title = "Volcano Plot",
        subtitle = paste0(contrast_name, "\n", "FDR < ", fdr_threshold, " & |log2FC| > ", logfc_threshold),
        x = "log2 Fold Change",
        y = "-log10(FDR)",
        color = "Classification"
      ) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    
    if (nrow(top_genes) > 0) {
      p <- p + geom_text_repel(data = top_genes, aes(label = gene),
                              size = 3, max.overlaps = 20)
    }
    
    # 4. Volcano plot
    ggsave(file.path(contrast_dir, "volcano_plot.png"),
           plot = p, width = 10, height = 8, dpi = 300)
    cat("    ✓ Volcano plot saved: volcano_plot.png\n")
    
    # 5. Save contrast-specific summary
    contrast_summary <- data.frame(
      Contrast = contrast_name,
      Total_Genes = nrow(all_results),
      Significant_FDR05_FC1 = n_sig,
      Upregulated = n_up,
      Downregulated = n_down,
      FDR_Threshold = fdr_threshold,
      LogFC_Threshold = logfc_threshold,
      stringsAsFactors = FALSE
    )
    write.csv(contrast_summary,
              file.path(contrast_dir, "contrast_summary.csv"),
              row.names = FALSE, quote = FALSE)
    cat("    ✓ Contrast summary saved\n")
    
    # Add to summary
    de_summary <- rbind(de_summary, data.frame(
      Contrast = contrast_name,
      Total_Genes = nrow(all_results),
      Significant_FDR05_FC1 = n_sig,
      Upregulated = n_up,
      Downregulated = n_down,
      stringsAsFactors = FALSE
    ))
  }
}

# === SAVE GENERAL COMBINED RESULTS (MAIN LIMMA FOLDER) ===
cat("\n=== Saving General Combined Results to Main Limma Folder ===\n")

# Combined results from all contrasts
if (length(all_results_combined) > 0) {
  combined_df <- do.call(rbind, all_results_combined)
  write.csv(combined_df,
            file.path(args$output_dir, "all_contrasts_combined_results.csv"),
            row.names = FALSE, quote = FALSE)
  cat("✓ Combined results saved: all_contrasts_combined_results.csv\n")
  cat("  Total rows:", nrow(combined_df), "(all genes × all contrasts)\n")
} else {
  cat("ℹ No contrasts processed - no combined results to save\n")
}

# Overall analysis summary
write.csv(de_summary,
          file.path(args$output_dir, "analysis_summary.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Analysis summary saved: analysis_summary.csv\n")

cat("\n=== LIMMA Analysis Complete ===\n")
cat("📁 Main results folder:", args$output_dir, "\n")
cat("📊 General results: normalized counts, summary, combined results\n")
cat("📁 Contrast-specific subfolders: individual results per contrast\n\n")
cat("📋 Summary of all contrasts:\n")
if (nrow(de_summary) > 0) {
  print(de_summary)
} else {
  cat("  No contrasts were processed.\n")
}
cat("\n")

# === SAVE R OBJECTS FOR REPRODUCIBILITY (MAIN LIMMA FOLDER) ===
cat("💾 Saving R objects for reproducibility...\n")
saveRDS(fit2, file.path(args$output_dir, "limma_fit_object.rds"))
saveRDS(v, file.path(args$output_dir, "voom_object.rds"))
saveRDS(cm, file.path(args$output_dir, "contrast_matrix.rds"))
saveRDS(dge, file.path(args$output_dir, "dge_object.rds"))
cat("✓ R objects saved for future analysis and reproducibility\n\n")

cat("🎉 LIMMA-voom differential expression analysis completed successfully!\n")
cat("📂 Check the results in:", args$output_dir, "\n")