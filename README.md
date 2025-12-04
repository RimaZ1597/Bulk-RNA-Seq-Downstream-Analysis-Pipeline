# Bulk-RNA-Seq-Downstream-Analysis-Pipeline
This pipeline performs reproducible RNA-Seq downstream analysis using DESeq2 with optional batch correction (ComBat-Seq), power analysis, differential expression, and pathway enrichment — all in a modular, containerised Nextflow workflow.

____

# Features

📦 Containerised with Docker/Singularity

📁 Configurable via analysis_config.yaml or interactive interactive_setup.py

🔄 Supports batch effect correction with ComBat-Seq

📊 Differential expression via DESeq2, limma-voom, or edgeR

🧬 Pathway enrichment (GSEA, clusterProfiler, KEGG, REACTOME)

📈 Power analysis with RNASeqPower

🧪 All results saved as unfiltered CSVs (no p-value/log2FC threshold applied)

✅ Reproducible, modular (DSL2) design


<img width="939" height="450" alt="workflow" src="https://github.com/user-attachments/assets/27122fcd-363c-40b9-a1e0-47ea18e63cc1" />

# Run Pipeline Using Interactive Setup

python interactive_setup.py
nextflow run main.nf -c nextflow.config --config_file config_summary.txt
