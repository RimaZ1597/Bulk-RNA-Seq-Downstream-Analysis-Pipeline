#!/bin/bash

# Test the fixed workflow pipeline using interactive setup configuration
# Usage: bash run_pipeline_test.sh

echo "🧬 Starting Fixed Workflow Pipeline Test"
echo "📁 Using configuration from: interactive_setup/analysis_config.yaml"
echo "🎯 Experiment: REPSOX Data Analysis"
echo ""

# Run the pipeline with interactive setup configuration
nextflow run main.nf \
    -profile standard \
    --config_file interactive_setup/analysis_config.yaml \
    --experiment_id "REPSOX_FixedWorkflow_Test" \
    -resume

echo ""
echo "✅ Pipeline test completed!"
echo "📊 Check results in: results/REPSOX_FixedWorkflow_Test/Results/"