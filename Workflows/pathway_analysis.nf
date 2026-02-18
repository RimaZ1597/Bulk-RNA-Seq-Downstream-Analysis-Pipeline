#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    PATHWAY ANALYSIS WORKFLOW
========================================================================================
    Performs FGSEA pathway analysis on limma differential expression results
    
    Input:  Limma results (per contrast) + subsystem genes + config
    Output: Pathway analysis results (HumanGEM + MSigDB) per contrast + combined summary
    
    Results Structure: {Experiment_ID}/Results/pathway_analysis/
      ├── all_contrasts_pathway_summary.csv
      └── {Experiment_ID}/Results/Limma/{Contrast}/pathway_analysis/
          ├── ranked_gene_list.csv
          ├── humangem_subsystem_pathways.csv
          ├── msigdb_*_pathways.csv
          ├── msigdb_all_collections_combined.csv
          └── pathway_analysis_summary.csv
========================================================================================
*/

include { PATHWAY_ANALYSIS; COMBINE_PATHWAY_RESULTS } from '../modules/local/pathway_analysis/main'

workflow PATHWAY_WORKFLOW {
    take:
    limma_results_ch  // [project_id, contrast_name, safe_contrast_name, limma_results.csv, config.yaml, subsystem_genes.csv]
    
    main:
    
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║            FGSEA Pathway Analysis Workflow                ║
    ╚═══════════════════════════════════════════════════════════╝
    """
    
    // Prepare input for pathway analysis (subsystem genes already included in tuples)
    pathway_input_ch = limma_results_ch
        .map { project_id, contrast_name, safe_contrast_name, limma_results, config_yaml, subsystem_file ->
            log.info "[PATHWAY_WORKFLOW] ▶ Processing pathway analysis for: ${contrast_name} (safe: ${safe_contrast_name})"
            tuple(project_id, contrast_name, safe_contrast_name, limma_results, config_yaml, subsystem_file)
        }
    
    // Add pathway analysis script
    pathway_script_ch = Channel.fromPath("${projectDir}/bin/pathway_analysis.R", checkIfExists: true)
    
    // Run pathway analysis for each contrast
    PATHWAY_ANALYSIS(pathway_input_ch, pathway_script_ch)
    
    // Collect all summaries for combining
    all_summaries_ch = PATHWAY_ANALYSIS.out.summary
        .map { project_id, contrast_name, summary_file -> 
            tuple(project_id, summary_file) 
        }
        .groupTuple()
    
    // Combine all pathway summaries
    COMBINE_PATHWAY_RESULTS(all_summaries_ch)
    
    emit:
    ranked_genes = PATHWAY_ANALYSIS.out.ranked_genes
    humangem_all = PATHWAY_ANALYSIS.out.humangem_all
    humangem_sig = PATHWAY_ANALYSIS.out.humangem_sig
    msigdb_collections = PATHWAY_ANALYSIS.out.msigdb_collections
    msigdb_combined = PATHWAY_ANALYSIS.out.msigdb_combined
    msigdb_sig = PATHWAY_ANALYSIS.out.msigdb_sig
    summary = PATHWAY_ANALYSIS.out.summary
    combined_summary = COMBINE_PATHWAY_RESULTS.out.combined_summary
    log_files = PATHWAY_ANALYSIS.out.log
}

workflow {
    // Standalone execution mode
    if (!params.project_id) error "ERROR: --project_id required"
    if (!params.limma_dir) error "ERROR: --limma_dir required (path to limma results directory)"
    if (!params.subsystem_genes) error "ERROR: --subsystem_genes required"
    if (!params.config) error "ERROR: --config required"
    
    // Parse limma directory to find all contrast results
    limma_dir = file(params.limma_dir)
    
    // Find all contrast directories
    contrast_dirs = Channel.fromPath("${limma_dir}/*", type: 'dir')
        .filter { it.name != 'limma' && !it.name.endsWith('.csv') && !it.name.endsWith('.rds') }
        .map { contrast_dir ->
            def contrast_name = contrast_dir.name
            def results_file = file("${contrast_dir}/all_genes_results.csv")
            
            if (!results_file.exists()) {
                log.warn "⚠ No all_genes_results.csv found in ${contrast_name}, skipping"
                return null
            }
            
            tuple(
                params.project_id,
                contrast_name,
                contrast_name,  // safe_contrast_name (already safe from limma)
                results_file,
                file(params.config)
            )
        }
        .filter { it != null }
    
    PATHWAY_WORKFLOW(
        contrast_dirs,
        file(params.subsystem_genes)
    )
}