#!/usr/bin/env Rscript

# pathway_analysis.R
# FGSEA Pathway Analysis for bulk RNA-seq
# Performs: HumanGEM subsystem analysis + MSigDB pathway analysis
# Output: CSV files only (contrast-specific + summary)

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(msigdbr)
  library(yaml)
  library(argparse)
})

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

parser <- ArgumentParser(description='FGSEA pathway analysis for RNA-seq')
parser$add_argument('--limma_results', required=TRUE,
                    help='Path to limma all_genes_results.csv for specific contrast')
parser$add_argument('--contrast_name', required=TRUE,
                    help='Contrast name (e.g., donor_19503_shear_stress_High_vs_Static)')
parser$add_argument('--subsystem_genes', required=TRUE,
                    help='Path to subsystem_genes.csv (HumanGEM)')
parser$add_argument('--config', required=TRUE,
                    help='Path to analysis_config.yaml')
parser$add_argument('--output_dir', default='.',
                    help='Output directory for results')
parser$add_argument('--species', default='Homo sapiens',
                    help='Species name for MSigDB')
parser$add_argument('--fdr_threshold', type='double', default=0.05,
                    help='FDR threshold for significant pathways')
parser$add_argument('--min_size', type='integer', default=15,
                    help='Minimum pathway size')
parser$add_argument('--max_size', type='integer', default=500,
                    help='Maximum pathway size')

args <- parser$parse_args()

cat("\n════════════════════════════════════════════════════════════\n")
cat("         FGSEA PATHWAY ANALYSIS PIPELINE\n")
cat("════════════════════════════════════════════════════════════\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Input Parameters:\n")
cat("  Limma results:", args$limma_results, "\n")
cat("  Contrast:", args$contrast_name, "\n")
cat("  Subsystem genes:", args$subsystem_genes, "\n")
cat("  Configuration:", args$config, "\n")
cat("  Output directory:", args$output_dir, "\n")
cat("  Species:", args$species, "\n")
cat("  FDR threshold:", args$fdr_threshold, "\n")
cat("  Pathway size:", args$min_size, "-", args$max_size, "\n\n")

# ============================================================================
# LOAD CONFIGURATION
# ============================================================================

cat("Loading analysis configuration...\n")
tryCatch({
  config <- yaml::read_yaml(args$config)
}, error = function(e) {
  stop("Failed to load configuration file: ", e$message)
})

project_name <- config$project_name
cat("  Project:", project_name, "\n\n")

# ============================================================================
# LOAD LIMMA RESULTS AND CREATE RANKED GENE LIST
# ============================================================================

cat("════════════════════════════════════════════════════════════\n")
cat("STEP 1: Loading Limma Results and Creating Ranked Gene List\n")
cat("════════════════════════════════════════════════════════════\n")

# Load Limma results
tryCatch({
  limma_df <- read.csv(args$limma_results, stringsAsFactors = FALSE)
  cat("✓ Loaded limma results:", nrow(limma_df), "genes\n")
}, error = function(e) {
  stop("Failed to load limma results: ", e$message)
})

# Validate required columns
required_cols <- c("gene", "logFC", "t", "P.Value", "adj.P.Val")
missing_cols <- setdiff(required_cols, colnames(limma_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in limma results: ", paste(missing_cols, collapse=", "))
}

# Create ranking metric: sign(logFC) * -log10(P.Value)
# This prioritizes genes with large fold changes and small p-values
limma_df$ranking_metric <- sign(limma_df$logFC) * -log10(limma_df$P.Value)

# Remove genes with NA or infinite values
limma_df <- limma_df[!is.na(limma_df$ranking_metric) & 
                     !is.infinite(limma_df$ranking_metric), ]

cat("  Genes after filtering:", nrow(limma_df), "\n")

# Create named vector for fgsea (gene -> ranking metric)
gene_ranks <- setNames(limma_df$ranking_metric, limma_df$gene)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

cat("  Ranked gene list created:", length(gene_ranks), "genes\n")
cat("    Top ranked gene:", names(gene_ranks)[1], "(", round(gene_ranks[1], 2), ")\n")
cat("    Bottom ranked gene:", names(gene_ranks)[length(gene_ranks)], 
    "(", round(gene_ranks[length(gene_ranks)], 2), ")\n\n")

# Save ranked gene list
ranked_genes_df <- data.frame(
  gene = names(gene_ranks),
  ranking_metric = gene_ranks,
  logFC = limma_df$logFC[match(names(gene_ranks), limma_df$gene)],
  P.Value = limma_df$P.Value[match(names(gene_ranks), limma_df$gene)],
  adj.P.Val = limma_df$adj.P.Val[match(names(gene_ranks), limma_df$gene)],
  stringsAsFactors = FALSE
)
write.csv(ranked_genes_df, 
          file.path(args$output_dir, "ranked_gene_list.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Saved: ranked_gene_list.csv\n\n")

# ============================================================================
# HUMANGEM SUBSYSTEM PATHWAY ANALYSIS
# ============================================================================

cat("════════════════════════════════════════════════════════════\n")
cat("STEP 2: HumanGEM Subsystem Pathway Analysis\n")
cat("════════════════════════════════════════════════════════════\n")

# Load subsystem genes
tryCatch({
  subsystem_df <- read.csv(args$subsystem_genes, stringsAsFactors = FALSE)
  cat("✓ Loaded subsystem genes file\n")
  cat("  Columns:", paste(colnames(subsystem_df), collapse=", "), "\n")
}, error = function(e) {
  stop("Failed to load subsystem_genes.csv: ", e$message)
})

# Handle different input formats
# Format 1: subsystem, gene_symbol (individual genes)  
# Format 2: gs_name, gene_list (comma-separated genes)
if ("gs_name" %in% colnames(subsystem_df) && "gene_list" %in% colnames(subsystem_df)) {
  cat("  Format: gs_name + gene_list (comma-separated)\n")
  # Convert comma-separated format to individual genes
  subsystem_pathways <- list()
  for (i in 1:nrow(subsystem_df)) {
    pathway_name <- subsystem_df$gs_name[i]
    genes <- trimws(strsplit(subsystem_df$gene_list[i], ",")[[1]])
    genes <- genes[genes != "" & !is.na(genes)]
    if (length(genes) > 0) {
      subsystem_pathways[[pathway_name]] <- unique(genes)
    }
  }
} else {
  # Format 1: individual gene rows
  if (!"subsystem" %in% colnames(subsystem_df)) {
    if ("Subsystem" %in% colnames(subsystem_df)) {
      colnames(subsystem_df)[colnames(subsystem_df) == "Subsystem"] <- "subsystem"
    } else {
      stop("subsystem_genes.csv must have 'subsystem'/'Subsystem' or 'gs_name' column")
    }
  }

  if (!"gene_symbol" %in% colnames(subsystem_df)) {
    if ("gene" %in% colnames(subsystem_df)) {
      colnames(subsystem_df)[colnames(subsystem_df) == "gene"] <- "gene_symbol"
    } else if ("Gene" %in% colnames(subsystem_df)) {
      colnames(subsystem_df)[colnames(subsystem_df) == "Gene"] <- "gene_symbol"
    } else {
      stop("subsystem_genes.csv must have 'gene_symbol'/'gene'/'Gene' or 'gene_list' column")
    }
  }

  # Convert to gene set list
  subsystem_pathways <- split(subsystem_df$gene_symbol, subsystem_df$subsystem)
  subsystem_pathways <- lapply(subsystem_pathways, unique)
}

cat("  Total subsystems:", length(subsystem_pathways), "\n")
cat("  Subsystem sizes:\n")
pathway_sizes <- sapply(subsystem_pathways, length)
cat("    Min:", min(pathway_sizes), "\n")
cat("    Max:", max(pathway_sizes), "\n")
cat("    Median:", median(pathway_sizes), "\n\n")

# Run FGSEA on HumanGEM subsystems
cat("Running FGSEA on HumanGEM subsystems...\n")
set.seed(42)
tryCatch({
  subsystem_fgsea <- fgsea(
    pathways = subsystem_pathways,
    stats = gene_ranks,
    minSize = args$min_size,
    maxSize = args$max_size,
    nperm = 10000
  )
  
  cat("✓ FGSEA completed\n")
  cat("  Total pathways tested:", nrow(subsystem_fgsea), "\n")
  cat("  Significant pathways (FDR < ", args$fdr_threshold, "): ",
      sum(subsystem_fgsea$padj < args$fdr_threshold, na.rm = TRUE), "\n\n")
  
}, error = function(e) {
  stop("FGSEA failed on HumanGEM subsystems: ", e$message)
})

# Format results
subsystem_results <- subsystem_fgsea %>%
  as.data.frame() %>%
  arrange(padj, pval) %>%
  mutate(
    pathway = as.character(pathway),
    leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=", "))
  ) %>%
  select(pathway, pval, padj, ES, NES, size, leadingEdge)

# Save results
write.csv(subsystem_results,
          file.path(args$output_dir, "humangem_subsystem_pathways.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Saved: humangem_subsystem_pathways.csv\n")

# Save significant pathways only
sig_subsystem <- subsystem_results %>%
  filter(padj < args$fdr_threshold)

if (nrow(sig_subsystem) > 0) {
  write.csv(sig_subsystem,
            file.path(args$output_dir, "humangem_subsystem_pathways_significant.csv"),
            row.names = FALSE, quote = FALSE)
  cat("✓ Saved: humangem_subsystem_pathways_significant.csv\n")
  cat("  Significant subsystems:", nrow(sig_subsystem), "\n\n")
} else {
  cat("ℹ  No significant HumanGEM subsystems at FDR <", args$fdr_threshold, "\n\n")
}

# ============================================================================
# MSIGDB PATHWAY ANALYSIS
# ============================================================================

cat("════════════════════════════════════════════════════════════\n")
cat("STEP 3: MSigDB Pathway Analysis\n")
cat("════════════════════════════════════════════════════════════\n")

# MSigDB collections to analyze
msigdb_collections <- c(
  "H" = "Hallmark",
  "C2" = "Curated (KEGG, Reactome, BioCarta)",
  "C5" = "Gene Ontology",
  "C6" = "Oncogenic Signatures"
)

all_msigdb_results <- list()

for (collection_code in names(msigdb_collections)) {
  collection_name <- msigdb_collections[collection_code]
  
  cat("\n─────────────────────────────────────────────────────────────\n")
  cat("Processing:", collection_name, "(", collection_code, ")\n")
  cat("─────────────────────────────────────────────────────────────\n")
  
  # Load MSigDB gene sets
  tryCatch({
    msigdb_sets <- msigdbr(species = args$species, category = collection_code)
    cat("  Loaded", nrow(msigdb_sets), "gene set entries\n")
    
    # Convert to list format for fgsea
    msigdb_pathways <- split(msigdb_sets$gene_symbol, msigdb_sets$gs_name)
    msigdb_pathways <- lapply(msigdb_pathways, unique)
    
    cat("  Total gene sets:", length(msigdb_pathways), "\n")
    
    # Run FGSEA
    set.seed(42)
    msigdb_fgsea <- fgsea(
      pathways = msigdb_pathways,
      stats = gene_ranks,
      minSize = args$min_size,
      maxSize = args$max_size,
      nperm = 10000
    )
    
    # Format results
    msigdb_results <- msigdb_fgsea %>%
      as.data.frame() %>%
      arrange(padj, pval) %>%
      mutate(
        pathway = as.character(pathway),
        collection = collection_name,
        leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=", "))
      ) %>%
      select(pathway, collection, pval, padj, ES, NES, size, leadingEdge)
    
    # Save collection-specific results
    clean_name <- gsub("[^A-Za-z0-9]", "_", collection_name)
    write.csv(msigdb_results,
              file.path(args$output_dir, paste0("msigdb_", clean_name, "_pathways.csv")),
              row.names = FALSE, quote = FALSE)
    
    cat("  ✓ Saved: msigdb_", clean_name, "_pathways.csv\n")
    cat("  Total pathways:", nrow(msigdb_results), "\n")
    cat("  Significant (FDR < ", args$fdr_threshold, "): ",
        sum(msigdb_results$padj < args$fdr_threshold, na.rm = TRUE), "\n")
    
    # Save significant pathways
    sig_msigdb <- msigdb_results %>% filter(padj < args$fdr_threshold)
    if (nrow(sig_msigdb) > 0) {
      write.csv(sig_msigdb,
                file.path(args$output_dir, paste0("msigdb_", clean_name, 
"_pathways_significant.csv")),
                row.names = FALSE, quote = FALSE)
      cat("  ✓ Saved significant pathways:", nrow(sig_msigdb), "\n")
    }
    
    # Store for combined results
    all_msigdb_results[[collection_name]] <- msigdb_results
    
  }, error = function(e) {
    cat("  ⚠ Warning: Failed to process", collection_name, "-", e$message, "\n")
  })
}

# ============================================================================
# COMBINED SUMMARY
# ============================================================================

cat("\n════════════════════════════════════════════════════════════\n")
cat("STEP 4: Creating Combined Summary\n")
cat("════════════════════════════════════════════════════════════\n")

# Combine all MSigDB results
if (length(all_msigdb_results) > 0) {
  combined_msigdb <- do.call(rbind, all_msigdb_results)
  write.csv(combined_msigdb,
            file.path(args$output_dir, "msigdb_all_collections_combined.csv"),
            row.names = FALSE, quote = FALSE)
  cat("✓ Saved: msigdb_all_collections_combined.csv\n")
  
  # Significant pathways across all collections
  sig_combined <- combined_msigdb %>% 
    filter(padj < args$fdr_threshold) %>%
    arrange(padj)
  
  if (nrow(sig_combined) > 0) {
    write.csv(sig_combined,
              file.path(args$output_dir, "msigdb_all_collections_significant.csv"),
              row.names = FALSE, quote = FALSE)
    cat("✓ Saved: msigdb_all_collections_significant.csv\n")
    cat("  Total significant pathways:", nrow(sig_combined), "\n")
  }
}

# Create pathway analysis summary
summary_data <- data.frame(
  analysis_type = c("HumanGEM Subsystems", names(msigdb_collections)),
  total_pathways = c(
    nrow(subsystem_results),
    sapply(all_msigdb_results, nrow)
  ),
  significant_pathways = c(
    sum(subsystem_results$padj < args$fdr_threshold, na.rm = TRUE),
    sapply(all_msigdb_results, function(x) sum(x$padj < args$fdr_threshold, na.rm = TRUE))
  ),
  fdr_threshold = args$fdr_threshold,
  min_pathway_size = args$min_size,
  max_pathway_size = args$max_size,
  stringsAsFactors = FALSE
)

write.csv(summary_data,
          file.path(args$output_dir, "pathway_analysis_summary.csv"),
          row.names = FALSE, quote = FALSE)
cat("✓ Saved: pathway_analysis_summary.csv\n\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("════════════════════════════════════════════════════════════\n")
cat("PATHWAY ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════\n")
cat("Contrast:", args$contrast_name, "\n")
cat("Output directory:", args$output_dir, "\n\n")

cat("Results Summary:\n")
cat("  HumanGEM Subsystems:\n")
cat("    Total pathways:", nrow(subsystem_results), "\n")
cat("    Significant:", sum(subsystem_results$padj < args$fdr_threshold, na.rm = TRUE), "\n\n")

for (collection_name in names(all_msigdb_results)) {
  results <- all_msigdb_results[[collection_name]]
  cat("  ", collection_name, ":\n")
  cat("    Total pathways:", nrow(results), "\n")
  cat("    Significant:", sum(results$padj < args$fdr_threshold, na.rm = TRUE), "\n\n")
}

cat("Output Files:\n")
cat("  ✓ ranked_gene_list.csv\n")
cat("  ✓ humangem_subsystem_pathways.csv\n")
cat("  ✓ msigdb_*_pathways.csv (per collection)\n")
cat("  ✓ msigdb_all_collections_combined.csv\n")
cat("  ✓ pathway_analysis_summary.csv\n")
cat("════════════════════════════════════════════════════════════\n\n")

cat("🎉 Pathway analysis completed successfully!\n")