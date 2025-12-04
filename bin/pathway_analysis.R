#!/usr/bin/env Rscript

# Enhanced Pathway Analysis - Separate CSV files for each collection (Multi-Dataset Support)
# Usage: Rscript pathway_analysis.R --deseq2_results <dir> --config <file> --output_dir <dir>

suppressPackageStartupMessages({
    library(fgsea)           # Fast GSEA
    library(msigdbr)         # MSigDB gene sets
    library(dplyr)           # Data manipulation
    library(yaml)            # Configuration
    library(argparse)        # Command line arguments
})

# Parse command line arguments
parser <- ArgumentParser(description='Enhanced pathway analysis with separate collection outputs and multi-dataset support')
parser$add_argument('--deseq2_results', required=TRUE, help='DESeq2 results directory')
parser$add_argument('--config', required=TRUE, help='Configuration YAML file')
parser$add_argument('--output_dir', default='.', help='Output directory')
args <- parser$parse_args()

cat("=== Enhanced Pathway Analysis Pipeline with Separate Collection Outputs ===\n")
cat("Parameters:\n")
cat("  DESeq2 results directory:", args$deseq2_results, "\n")
cat("  Config file:", args$config, "\n")
cat("  Output directory:", args$output_dir, "\n\n")

# Load configuration
tryCatch({
    config <- yaml::read_yaml(args$config)
    cat("✓ Configuration loaded successfully\n")
}, error = function(e) {
    stop("Failed to load configuration: ", e$message)
})

#####################################
## Configuration Defaults
#####################################

# Set defaults with config override
SPECIES <- config$species %||% "Homo sapiens"
MIN_PATHWAY_SIZE <- config$min_pathway_size %||% 10
MAX_PATHWAY_SIZE <- config$max_pathway_size %||% 1000
P_VALUE_THRESHOLD <- config$p_value_threshold %||% 1.0  # No filtering by default

cat("Configuration:\n")
cat("  Species:", SPECIES, "\n")
cat("  Min pathway size:", MIN_PATHWAY_SIZE, "\n")
cat("  Max pathway size:", MAX_PATHWAY_SIZE, "\n")
cat("  P-value threshold:", P_VALUE_THRESHOLD, "\n\n")

#####################################
## Helper Functions
#####################################

validate_deseq_results <- function(deseq_results, contrast_name) {
    cat("Validating DESeq2 results...\n")
    # Check for required columns
    required_cols <- c("stat")
    missing_cols <- setdiff(required_cols, colnames(deseq_results))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns in", contrast_name, ":", paste(missing_cols, collapse = ", ")))
    }
    # Check for valid statistics
    valid_stats <- sum(!is.na(deseq_results$stat))
    if (valid_stats < 50) {
        warning(paste("Very few genes with valid statistics in", contrast_name, ":", valid_stats))
    }
    cat(sprintf("  Valid genes with statistics: %d\n", valid_stats))
    return(TRUE)
}

detect_gene_symbols <- function(deseq_results) {
    cat("Detecting gene symbols...\n")
    # Common gene identifier column names
    gene_cols <- c("gene", "Gene", "gene_symbol", "symbol", "SYMBOL",
                   "GeneName", "gene_name", "hgnc_symbol", "external_gene_name")
    # Check for gene column in data
    gene_col <- intersect(gene_cols, colnames(deseq_results))[1]
    if (!is.na(gene_col)) {
        cat(sprintf("  Found gene symbols in column: %s\n", gene_col))
        rownames(deseq_results) <- make.unique(as.character(deseq_results[[gene_col]]))
        # Remove the gene column to avoid duplication
        deseq_results[[gene_col]] <- NULL
        return(deseq_results)
    }
    # Check if rownames look like gene symbols
    sample_names <- rownames(deseq_results)[1:min(10, nrow(deseq_results))]
    # Pattern for human gene symbols (start with letter, contain letters/numbers/hyphens)
    if (all(grepl("^[A-Za-z][A-Za-z0-9.-]*$", sample_names)) &&
        !all(grepl("^ENSG", sample_names))) {  # Not Ensembl IDs
        cat("  Using existing rownames as gene symbols\n")
        return(deseq_results)
    }
    # Check for Ensembl IDs in rownames
    if (all(grepl("^ENS", sample_names))) {
        warning("Data appears to contain Ensembl IDs. Gene symbol conversion may be needed.")
        cat("  Warning: Ensembl IDs detected - pathway analysis may have reduced coverage\n")
        return(deseq_results)
    }
    stop("Cannot identify valid gene symbols in the data. Please ensure gene symbols are available.")
}

find_deseq_file <- function(contrast_dir) {
    # Possible file locations and patterns
    possible_locations <- c(
        file.path(contrast_dir, "Differential_Expression_Analysis", "DESeq2_results.csv"),
        file.path(contrast_dir, "DESeq2_results.csv"),
        file.path(contrast_dir, "deseq2_results.csv"),
        file.path(contrast_dir, "results.csv")
    )
    # Check exact matches first
    for (file_path in possible_locations) {
        if (file.exists(file_path)) {
            cat(sprintf("  Found DESeq2 file: %s\n", basename(file_path)))
            return(file_path)
        }
    }
    # Search for pattern matches
    pattern_files <- list.files(contrast_dir,
                               pattern = "(deseq|results|differential).*\\.csv$",
                               recursive = TRUE,
                               full.names = TRUE,
                               ignore.case = TRUE)
    if (length(pattern_files) > 0) {
        cat(sprintf("  Found DESeq2 file by pattern: %s\n", basename(pattern_files[1])))
        return(pattern_files[1])
    }
    return(NULL)
}

get_msigdb_pathways <- function(species = "Homo sapiens") {
    cat(sprintf("Loading MSigDB pathway collections for %s...\n", species))
    
    pathway_collections <- NULL  # Initialize to avoid "not found" error
    
    tryCatch({
        # Try with new API (collection/subcollection)
        tryCatch({
            # Hallmark gene sets
            pathways_h <- msigdbr(species = species, collection = "H") %>%
                split(x = .$gene_symbol, f = .$gs_name)
            
            # C2 - Canonical Pathways (General)
            pathways_c2_cp <- msigdbr(species = species, collection = "C2", subcollection = "CP") %>%
                split(x = .$gene_symbol, f = .$gs_name)
            
            # C2 - KEGG Pathways
            pathways_c2_kegg <- tryCatch({
                msigdbr(species = species, collection = "C2", subcollection = "CP:KEGG") %>%
                    split(x = .$gene_symbol, f = .$gs_name)
            }, error = function(e) {
                cat("KEGG subcollection failed, using C2 canonical pathways\n")
                pathways_c2_cp  # Use canonical pathways as fallback
            })
            
            # C2 - Reactome Pathways
            pathways_c2_reactome <- tryCatch({
                msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME") %>%
                    split(x = .$gene_symbol, f = .$gs_name)
            }, error = function(e) {
                cat("Reactome subcollection failed\n")
                list()
            })
            
            # C5 - GO Gene Sets
            pathways_c5_bp <- tryCatch({
                msigdbr(species = species, collection = "C5", subcollection = "GO:BP") %>%
                    split(x = .$gene_symbol, f = .$gs_name)
            }, error = function(e) {
                cat("GO:BP subcollection failed\n")
                list()
            })
            
            pathway_collections <- list(
                "Hallmark" = pathways_h,
                "C2_Canonical_Pathways" = pathways_c2_cp,
                "C2_KEGG" = pathways_c2_kegg,
                "C2_Reactome" = pathways_c2_reactome,
                "C5_GO_Biological_Process" = pathways_c5_bp
            )
            
        }, error = function(e) {
            cat("New API failed, trying old API...\n")
            
            # Fallback to old API if new one fails
            tryCatch({
                pathways_h <- msigdbr(species = species, category = "H") %>%
                    split(x = .$gene_symbol, f = .$gs_name)
                
                pathways_c2 <- msigdbr(species = species, category = "C2") %>%
                    split(x = .$gene_symbol, f = .$gs_name)
                
                pathways_c5 <- msigdbr(species = species, category = "C5") %>%
                    split(x = .$gene_symbol, f = .$gs_name)
                
                pathway_collections <- list(
                    "Hallmark" = pathways_h,
                    "C2_All_Pathways" = pathways_c2,
                    "C5_Gene_Ontology" = pathways_c5
                )
                
            }, error = function(e2) {
                cat("Both APIs failed, using minimal fallback...\n")
                
                # Ultimate fallback - just Hallmark
                pathways_h <- tryCatch({
                    if ("collection" %in% names(formals(msigdbr))) {
                        msigdbr(species = species, collection = "H")
                    } else {
                        msigdbr(species = species, category = "H")
                    }
                }, error = function(e3) {
                    msigdbr(species = "Homo sapiens", category = "H")  # Force human as last resort
                }) %>%
                    split(x = .$gene_symbol, f = .$gs_name)
                
                pathway_collections <- list(
                    "Hallmark" = pathways_h
                )
            })
        })
        
    }, error = function(e) {
        cat("All pathway loading attempts failed\n")
        cat(sprintf("Error details: %s\n", e$message))
        
        # Create empty collections to prevent script failure
        pathway_collections <- list(
            "Hallmark" = list(),
            "C2_Pathways" = list()
        )
        
        warning("Could not load MSigDB pathways - continuing with empty collections")
    })
    
    # Ensure pathway_collections is not NULL
    if (is.null(pathway_collections)) {
        cat("pathway_collections is NULL, creating empty collections\n")
        pathway_collections <- list(
            "Hallmark" = list(),
            "C2_Pathways" = list()
        )
    }
    
    # Filter out empty collections
    pathway_collections <- pathway_collections[sapply(pathway_collections, length) > 0]
    
    cat("Pathway Collections Summary:\n")
    if (length(pathway_collections) > 0) {
        for (name in names(pathway_collections)) {
            cat(sprintf("  %s: %d pathways\n", name, length(pathway_collections[[name]])))
        }
    } else {
        cat("  No pathway collections loaded successfully\n")
        warning("No pathway collections available - pathway analysis may be limited")
    }
    cat("\n")
    
    return(pathway_collections)
}

prepare_ranks <- function(deseq_result) {
    cat("Preparing gene ranks...\n")
    # Use stat column for ranking (t-statistic or Wald statistic)
    ranks <- deseq_result$stat
    names(ranks) <- rownames(deseq_result)
    # Remove NA values
    ranks <- ranks[!is.na(ranks)]
    # Remove duplicated gene names (keep first occurrence)
    ranks <- ranks[!duplicated(names(ranks))]
    # Sort in decreasing order (most upregulated first)
    ranks <- sort(ranks, decreasing = TRUE)
    cat(sprintf("  Prepared ranks for %d genes\n", length(ranks)))
    cat(sprintf("  Range: %.3f to %.3f\n", max(ranks), min(ranks)))
    return(ranks)
}

# Enhanced function to add pathway metadata and enrichment information
enhance_pathway_results <- function(results_df, collection_name, contrast_name, ranks) {
    if (nrow(results_df) == 0) return(results_df)
    # Add basic metadata
    results_df <- results_df %>%
        mutate(
            collection = collection_name,
            contrast = contrast_name,
            abs_NES = abs(NES),
            direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
            significance_level = case_when(
                pval < 0.001 ~ "p < 0.001",
                pval < 0.01 ~ "p < 0.01",
                pval < 0.05 ~ "p < 0.05",
                TRUE ~ "NS"
            ),
            # Add FDR categories
            fdr_category = case_when(
                padj < 0.001 ~ "FDR < 0.001",
                padj < 0.01 ~ "FDR < 0.01",
                padj < 0.05 ~ "FDR < 0.05",
                padj < 0.1 ~ "FDR < 0.1",
                TRUE ~ "FDR >= 0.1"
            ),
            # Effect size categories
            effect_size_category = case_when(
                abs_NES >= 2.5 ~ "Very Large",
                abs_NES >= 2.0 ~ "Large",
                abs_NES >= 1.5 ~ "Medium",
                abs_NES >= 1.0 ~ "Small",
                TRUE ~ "Very Small"
            ),
            # Leading edge statistics
            leading_edge_fraction = sapply(1:nrow(results_df), function(i) {
                if (results_df$size[i] > 0) {
                    length(unlist(strsplit(results_df$leadingEdge[i], ";"))) / results_df$size[i]
                } else {
                    0
                }
            })
        )
    
    # Convert leadingEdge list to semicolon-separated string if it's still a list
    if (is.list(results_df$leadingEdge)) {
        results_df$leadingEdge <- sapply(results_df$leadingEdge, function(x) paste(x, collapse = ";"))
    }
    
    # Add pathway category based on pathway name patterns
    results_df <- results_df %>%
        mutate(
            pathway_category = case_when(
                # Hallmark categories
                grepl("INFLAMMATORY", pathway, ignore.case = TRUE) ~ "Inflammation",
                grepl("IMMUNE", pathway, ignore.case = TRUE) ~ "Immune Response",
                grepl("APOPTOSIS", pathway, ignore.case = TRUE) ~ "Cell Death",
                grepl("CELL_CYCLE|MITOTIC", pathway, ignore.case = TRUE) ~ "Cell Cycle",
                grepl("METABOLISM|METABOLIC", pathway, ignore.case = TRUE) ~ "Metabolism",
                grepl("SIGNALING", pathway, ignore.case = TRUE) ~ "Signaling",
                grepl("DEVELOPMENT", pathway, ignore.case = TRUE) ~ "Development",
                grepl("DNA_REPAIR", pathway, ignore.case = TRUE) ~ "DNA Repair",
                grepl("HYPOXIA", pathway, ignore.case = TRUE) ~ "Hypoxia Response",
                grepl("ANGIOGENESIS", pathway, ignore.case = TRUE) ~ "Angiogenesis",
                # GO categories
                grepl("GO_", pathway) & grepl("REGULATION", pathway, ignore.case = TRUE) ~ "Regulation",
                grepl("GO_", pathway) & grepl("TRANSPORT", pathway, ignore.case = TRUE) ~ "Transport",
                grepl("GO_", pathway) & grepl("BINDING", pathway, ignore.case = TRUE) ~ "Binding",
                grepl("GO_", pathway) & grepl("ACTIVITY", pathway, ignore.case = TRUE) ~ "Enzymatic Activity",
                # KEGG categories  
                grepl("KEGG_", pathway) & grepl("PATHWAY", pathway, ignore.case = TRUE) ~ "KEGG Pathway",
                # Reactome categories
                grepl("REACTOME_", pathway) ~ "Reactome Pathway",
                # C3 categories
                grepl("MIR_", pathway) ~ "MicroRNA Target",
                grepl("_TF_", pathway) ~ "Transcription Factor Target",
                # Default
                TRUE ~ "Other"
            )
        ) %>%
        # Reorder columns for better readability
        select(pathway, collection, contrast, pval, padj, NES, abs_NES, direction,
               significance_level, fdr_category, effect_size_category,
               size, leading_edge_fraction, pathway_category, leadingEdge, everything())
    
    return(results_df)
}

run_pathway_analysis <- function(ranks_list, pathway_collections, contrast_name, output_dir) {
    cat("Running enhanced pathway analysis...\n")
    results_list <- list()
    collection_stats <- list()
    
    for (collection in names(pathway_collections)) {
        cat(sprintf("  Processing collection: %s\n", collection))
        tryCatch({
            # Run fgsea
            results <- fgsea(
                pathways = pathway_collections[[collection]],
                stats = ranks_list,
                minSize = MIN_PATHWAY_SIZE,
                maxSize = MAX_PATHWAY_SIZE,
                eps = 0.0,
                nproc = 1  # Single thread for stability
            )
            
            if (nrow(results) > 0) {
                # Enhance results with additional information
                results_df <- as.data.frame(results)
                results_enhanced <- enhance_pathway_results(results_df, collection, contrast_name, ranks_list)
                # Sort by p-value, then by absolute NES
                results_enhanced <- results_enhanced %>%
                    arrange(pval, desc(abs_NES))
                results_list[[collection]] <- results_enhanced
                
                # Collection statistics
                collection_stats[[collection]] <- list(
                    total_pathways = nrow(results_enhanced),
                    significant_005 = sum(results_enhanced$pval < 0.05, na.rm = TRUE),
                    significant_001 = sum(results_enhanced$pval < 0.01, na.rm = TRUE),
                    upregulated = sum(results_enhanced$NES > 0, na.rm = TRUE),
                    downregulated = sum(results_enhanced$NES < 0, na.rm = TRUE),
                    large_effect = sum(results_enhanced$abs_NES >= 2.0, na.rm = TRUE)
                )
                
                cat(sprintf("    Found %d pathways (%d significant at p<0.05)\n",
                           nrow(results_enhanced),
                           sum(results_enhanced$pval < 0.05, na.rm = TRUE)))
            } else {
                cat(sprintf("    No pathways found for %s\n", collection))
            }
        }, error = function(e) {
            cat(sprintf("    Error in processing collection %s: %s\n", collection, e$message))
        })
        # Memory management
        gc()
    }
    
    if (length(results_list) == 0) {
        warning(paste("No pathway results generated for", contrast_name))
        return(list(combined = data.frame(), individual = list(), stats = list()))
    }
    
    # Save individual collection results with enhanced information
    for (collection in names(results_list)) {
        if (nrow(results_list[[collection]]) > 0) {
            collection_clean <- gsub("[: ()]", "_", collection)  # Clean filename
            
            # All results for this collection
            write.csv(
                results_list[[collection]],
                file = file.path(output_dir, paste0(collection_clean, "_all_pathways.csv")),
                row.names = FALSE
            )
            
            # Significant results (p < 0.05)
            sig_results <- results_list[[collection]] %>%
                filter(pval < 0.05) %>%
                arrange(pval, desc(abs_NES))
            
            if (nrow(sig_results) > 0) {
                write.csv(
                    sig_results,
                    file = file.path(output_dir, paste0(collection_clean, "_significant_pathways.csv")),
                    row.names = FALSE
                )
            }
            
            # Top pathways (top 50 by p-value)
            top_results <- results_list[[collection]] %>%
                arrange(pval) %>%
                head(50)
            
            write.csv(
                top_results,
                file = file.path(output_dir, paste0(collection_clean, "_top50_pathways.csv")),
                row.names = FALSE
            )
        }
    }
    
    # Combine all results
    combined_results <- bind_rows(results_list)
    
    # Save combined results
    write.csv(
        combined_results,
        file = file.path(output_dir, "all_collections_combined.csv"),
        row.names = FALSE
    )
    
    # Save significant combined results
    significant_combined <- combined_results %>%
        filter(pval < 0.05) %>%
        arrange(pval, desc(abs_NES))
    
    if (nrow(significant_combined) > 0) {
        write.csv(
            significant_combined,
            file = file.path(output_dir, "all_collections_significant.csv"),
            row.names = FALSE
        )
    }
    
    return(list(
        combined = combined_results,
        individual = results_list,
        stats = collection_stats
    ))
}

create_enhanced_summary <- function(pathway_results, contrast_name, output_dir) {
    cat("Creating enhanced summary statistics...\n")
    combined_results <- pathway_results$combined
    collection_stats <- pathway_results$stats
    
    if (nrow(combined_results) == 0) {
        cat("No results to summarize\n")
        return(NULL)
    }
    
    # Overall summary
    overall_summary <- data.frame(
        Contrast = contrast_name,
        Total_Pathways_Tested = nrow(combined_results),
        Significant_001 = sum(combined_results$pval < 0.001, na.rm = TRUE),
        Significant_01 = sum(combined_results$pval < 0.01, na.rm = TRUE),
        Significant_05 = sum(combined_results$pval < 0.05, na.rm = TRUE),
        FDR_001 = sum(combined_results$padj < 0.001, na.rm = TRUE),
        FDR_01 = sum(combined_results$padj < 0.01, na.rm = TRUE),
        FDR_05 = sum(combined_results$padj < 0.05, na.rm = TRUE),
        Upregulated_Pathways = sum(combined_results$NES > 0, na.rm = TRUE),
        Downregulated_Pathways = sum(combined_results$NES < 0, na.rm = TRUE),
        Large_Effect_Size = sum(combined_results$abs_NES >= 2.0, na.rm = TRUE),
        Medium_Effect_Size = sum(combined_results$abs_NES >= 1.5 & combined_results$abs_NES < 2.0, na.rm = TRUE),
        Analysis_Timestamp = as.character(Sys.time()),
        stringsAsFactors = FALSE
    )
    
    # Collection-wise summary
    collection_summary <- do.call(rbind, lapply(names(collection_stats), function(coll) {
        stats <- collection_stats[[coll]]
        data.frame(
            Contrast = contrast_name,
            Collection = coll,
            Total_Pathways = stats$total_pathways,
            Significant_05 = stats$significant_005,
            Significant_01 = stats$significant_001,
            Upregulated = stats$upregulated,
            Downregulated = stats$downregulated,
            Large_Effect_Size = stats$large_effect,
            Success_Rate = round(stats$significant_005 / stats$total_pathways * 100, 2),
            stringsAsFactors = FALSE
        )
    }))
    
    # Category summary (if pathway_category column exists)
    if ("pathway_category" %in% colnames(combined_results)) {
        category_summary <- combined_results %>%
            filter(pval < 0.05) %>%
            group_by(pathway_category, direction) %>%
            summarise(
                Count = n(),
                Mean_NES = mean(abs_NES, na.rm = TRUE),
                Mean_pval = mean(pval, na.rm = TRUE),
                .groups = 'drop'
            ) %>%
            mutate(Contrast = contrast_name) %>%
            arrange(desc(Count))
        
        write.csv(category_summary,
                 file.path(output_dir, "pathway_category_summary.csv"),
                 row.names = FALSE)
    }
    
    # Save summaries
    write.csv(overall_summary,
             file.path(output_dir, "overall_pathway_summary.csv"),
             row.names = FALSE)
    
    write.csv(collection_summary,
             file.path(output_dir, "collection_wise_summary.csv"),
             row.names = FALSE)
    
    return(list(
        overall = overall_summary,
        collections = collection_summary,
        categories = if(exists("category_summary")) category_summary else NULL
    ))
}

process_contrast <- function(contrast_dir, contrast_name, pathway_collections, config) {
    cat(sprintf("\n=== Processing Contrast: %s ===\n", contrast_name))
    
    # Find DESeq2 results file
    deseq_file <- find_deseq_file(contrast_dir)
    if (is.null(deseq_file)) {
        cat(sprintf("No DESeq2 results found for %s, skipping\n", contrast_name))
        return(NULL)
    }
    
    # Create output directory
    pathway_output_dir <- file.path(contrast_dir, "Pathway_Analysis")
    dir.create(pathway_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    tryCatch({
        # Load DESeq2 results
        deseq_results <- read.csv(deseq_file, row.names = NULL, stringsAsFactors = FALSE)
        cat(sprintf("Loaded data with %d rows and %d columns\n", nrow(deseq_results), ncol(deseq_results)))
        
        # Validate and process gene symbols
        validate_deseq_results(deseq_results, contrast_name)
        deseq_results <- detect_gene_symbols(deseq_results)
        
        # Prepare ranks
        ranks <- prepare_ranks(deseq_results)
        if (length(ranks) < 100) {
            warning(paste("Very few genes available for analysis in", contrast_name))
        }
        
        # Run enhanced pathway analysis
        pathway_results <- run_pathway_analysis(ranks, pathway_collections, contrast_name, pathway_output_dir)
        
        if (nrow(pathway_results$combined) > 0) {
            # Create enhanced summary statistics
            summary_stats <- create_enhanced_summary(pathway_results, contrast_name, pathway_output_dir)
            
            cat(sprintf("Enhanced pathway analysis completed for %s:\n", contrast_name))
            cat(sprintf("  Total pathways tested: %d\n", nrow(pathway_results$combined)))
            cat(sprintf("  Significant (p<0.05): %d\n", sum(pathway_results$combined$pval < 0.05, na.rm = TRUE)))
            cat(sprintf("  Significant (p<0.01): %d\n", sum(pathway_results$combined$pval < 0.01, na.rm = TRUE)))
            cat(sprintf("  FDR significant (q<0.05): %d\n", sum(pathway_results$combined$padj < 0.05, na.rm = TRUE)))
            cat(sprintf("  Large effect size (|NES|>=2): %d\n", sum(pathway_results$combined$abs_NES >= 2.0, na.rm = TRUE)))
        } else {
            cat(sprintf("No pathway results generated for %s\n", contrast_name))
        }
        
        return(pathway_results)
        
    }, error = function(e) {
        cat(sprintf("Error processing contrast %s: %s\n", contrast_name, e$message))
        cat(sprintf("Error details: %s\n", toString(e)))
        return(NULL)
    })
}

#####################################
## Main Analysis
#####################################

cat("\n=== Discovering Contrasts ===\n")

# Find contrast directories
contrast_base_dir <- file.path(args$deseq2_results, "03_Differential_Expression_&_Pathway_Analysis")

if (dir.exists(contrast_base_dir)) {
    cat("Using structured directory format\n")
    contrast_dirs <- list.dirs(contrast_base_dir, recursive = FALSE)
} else {
    cat("Using direct results directory format\n")
    contrast_dirs <- list.dirs(args$deseq2_results, recursive = FALSE)
    # If no subdirectories, try the root directory itself
    if (length(contrast_dirs) == 0) {
        contrast_dirs <- args$deseq2_results
    }
}

cat(sprintf("Found %d potential contrast directories\n", length(contrast_dirs)))

# Load pathway collections
pathway_collections <- get_msigdb_pathways(SPECIES)

# Process contrasts
processed_contrasts <- list()
failed_contrasts <- character()
successful_contrasts <- 0

# Process all found contrast directories
cat("Processing all found contrast directories\n")
for (contrast_dir in contrast_dirs) {
    contrast_name <- basename(contrast_dir)
    
    # Skip if directory name suggests it's not a contrast
    if (contrast_name %in% c("03_Differential_Expression_&_Pathway_Analysis", "results", "output")) {
        next
    }
    
    result <- process_contrast(contrast_dir, contrast_name, pathway_collections, config)
    if (!is.null(result)) {
        processed_contrasts[[contrast_name]] <- result
        successful_contrasts <- successful_contrasts + 1
        } else {
        failed_contrasts <- c(failed_contrasts, contrast_name)
    }
}

#####################################
## Create Global Summary
#####################################

cat("\n=== Creating Global Analysis Summary ===\n")

# Create global summary across all contrasts
if (length(processed_contrasts) > 0) {
    global_summary <- data.frame()
    global_collection_summary <- data.frame()
    
    for (contrast_name in names(processed_contrasts)) {
        contrast_results <- processed_contrasts[[contrast_name]]
        
        if (!is.null(contrast_results) && nrow(contrast_results$combined) > 0) {
            # Overall summary for this contrast
            contrast_summary <- data.frame(
                Contrast = contrast_name,
                Total_Pathways = nrow(contrast_results$combined),
                Significant_05 = sum(contrast_results$combined$pval < 0.05, na.rm = TRUE),
                Significant_01 = sum(contrast_results$combined$pval < 0.01, na.rm = TRUE),
                FDR_05 = sum(contrast_results$combined$padj < 0.05, na.rm = TRUE),
                FDR_01 = sum(contrast_results$combined$padj < 0.01, na.rm = TRUE),
                Upregulated = sum(contrast_results$combined$NES > 0, na.rm = TRUE),
                Downregulated = sum(contrast_results$combined$NES < 0, na.rm = TRUE),
                Large_Effect = sum(contrast_results$combined$abs_NES >= 2.0, na.rm = TRUE),
                Success_Rate = round(sum(contrast_results$combined$pval < 0.05, na.rm = TRUE) / nrow(contrast_results$combined) * 100, 2),
                stringsAsFactors = FALSE
            )
            global_summary <- rbind(global_summary, contrast_summary)
            
            # Collection-wise summary for this contrast
            for (collection_name in names(contrast_results$individual)) {
                collection_results <- contrast_results$individual[[collection_name]]
                if (!is.null(collection_results) && nrow(collection_results) > 0) {
                    collection_row <- data.frame(
                        Contrast = contrast_name,
                        Collection = collection_name,
                        Total_Pathways = nrow(collection_results),
                        Significant_05 = sum(collection_results$pval < 0.05, na.rm = TRUE),
                        Significant_01 = sum(collection_results$pval < 0.01, na.rm = TRUE),
                        FDR_05 = sum(collection_results$padj < 0.05, na.rm = TRUE),
                        Upregulated = sum(collection_results$NES > 0, na.rm = TRUE),
                        Downregulated = sum(collection_results$NES < 0, na.rm = TRUE),
                        Large_Effect = sum(collection_results$abs_NES >= 2.0, na.rm = TRUE),
                        stringsAsFactors = FALSE
                    )
                    global_collection_summary <- rbind(global_collection_summary, collection_row)
                }
            }
        }
    }
    
    # Save global summaries
    if (nrow(global_summary) > 0) {
        write.csv(global_summary,
                 file.path(args$output_dir, "global_pathway_summary.csv"),
                 row.names = FALSE)
        
        cat("Global Summary Statistics:\n")
        cat(sprintf("  Total contrasts processed: %d\n", nrow(global_summary)))
        cat(sprintf("  Average pathways per contrast: %.1f\n", mean(global_summary$Total_Pathways)))
        cat(sprintf("  Average significant pathways (p<0.05): %.1f\n", mean(global_summary$Significant_05)))
        cat(sprintf("  Average success rate: %.1f%%\n", mean(global_summary$Success_Rate)))
    }
    
    if (nrow(global_collection_summary) > 0) {
        write.csv(global_collection_summary,
                 file.path(args$output_dir, "global_collection_summary.csv"),
                 row.names = FALSE)
        
        # Collection performance summary
        collection_performance <- global_collection_summary %>%
            group_by(Collection) %>%
            summarise(
                Contrasts_Processed = n(),
                Avg_Total_Pathways = mean(Total_Pathways),
                Avg_Significant_05 = mean(Significant_05),
                Avg_FDR_05 = mean(FDR_05),
                Total_Significant = sum(Significant_05),
                .groups = 'drop'
            ) %>%
            arrange(desc(Total_Significant))
        
        write.csv(collection_performance,
                 file.path(args$output_dir, "collection_performance_summary.csv"),
                 row.names = FALSE)
        
        cat("\nTop performing collections:\n")
        for (i in 1:min(5, nrow(collection_performance))) {
            cat(sprintf("  %s: %d total significant pathways\n", 
                       collection_performance$Collection[i], 
                       collection_performance$Total_Significant[i]))
        }
    }
}

#####################################
## Final Report
#####################################

cat("\n=== Creating Final Analysis Report ===\n")

# Create comprehensive final report
final_report <- paste(
    "Enhanced Pathway Analysis Pipeline - Final Report",
    "===============================================",
    "",
    paste("Analysis completed:", Sys.time()),
    paste("Total contrasts found:", length(contrast_dirs)),
    paste("Successfully processed:", successful_contrasts),
    paste("Failed contrasts:", length(failed_contrasts)),
    "",
    "Configuration Used:",
    paste("- Species:", SPECIES),
    paste("- Min pathway size:", MIN_PATHWAY_SIZE),
    paste("- Max pathway size:", MAX_PATHWAY_SIZE),
    paste("- P-value threshold:", P_VALUE_THRESHOLD),
    "",
    "Pathway Collections Analyzed:",
    paste(sapply(names(pathway_collections), function(x) 
        paste("-", x, ":", length(pathway_collections[[x]]), "pathways")), 
        collapse = "\n"),
    "",
    if (length(failed_contrasts) > 0) {
        paste("Failed Contrasts:", paste(failed_contrasts, collapse = ", "))
    } else {
        "All contrasts processed successfully!"
    },
    "",
    "Output Structure:",
    "├── [Contrast_Name]/",
    "│   └── Pathway_Analysis/",
    "│       ├── [Collection]_all_pathways.csv",
    "│       ├── [Collection]_significant_pathways.csv",
    "│       ├── [Collection]_top50_pathways.csv",
    "│       ├── all_collections_combined.csv",
    "│       ├── all_collections_significant.csv",
    "│       ├── overall_pathway_summary.csv",
    "│       ├── collection_wise_summary.csv",
    "│       └── pathway_category_summary.csv",
    "",
    "Global Summary Files:",
    "├── global_pathway_summary.csv",
    "├── global_collection_summary.csv",
    "└── collection_performance_summary.csv",
    "",
    "Key Features:",
    "- Separate CSV files for each pathway collection",
    "- Enhanced pathway annotation with categories",
    "- Leading edge analysis and effect size classification",
    "- Multiple significance thresholds (p-value and FDR)",
    "- Comprehensive summary statistics",
    "- Global analysis across all contrasts",
    "",
    "Analysis completed successfully!",
    "",
    sep = "\n"
)

# Save final report
writeLines(final_report, file.path(args$output_dir, "pathway_analysis_final_report.txt"))

cat("\n" + "="*80)
cat("\n 🧬 ENHANCED PATHWAY ANALYSIS COMPLETED SUCCESSFULLY!")
cat("\n" + "="*80)
cat(sprintf("\n📊 Summary:\n"))
cat(sprintf("   ✅ Contrasts processed: %d/%d\n", successful_contrasts, length(contrast_dirs)))
cat(sprintf("   📁 Results organized by contrast and collection\n"))
cat(sprintf("   📋 Global summaries created\n"))

if (length(failed_contrasts) > 0) {
    cat(sprintf("   ⚠️  Failed contrasts: %d\n", length(failed_contrasts)))
    cat("   Failed contrast names:", paste(failed_contrasts, collapse = ", "), "\n")
}

cat(sprintf("\n📂 Results location: %s\n", args$output_dir))
cat(sprintf("📄 Final report: pathway_analysis_final_report.txt\n"))
cat(sprintf("🎯 Analysis features: Separate collection CSVs, enhanced annotations, global summaries\n"))

cat("\n" + "="*80 + "\n")

# Final memory cleanup
gc()

cat("Pipeline execution completed.\n")