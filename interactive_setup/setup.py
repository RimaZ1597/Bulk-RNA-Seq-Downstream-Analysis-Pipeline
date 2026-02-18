#!/usr/bin/env python3
"""
setup.py
========
Main orchestrator for downstream RNA-seq analysis configuration.
"""

import re
import yaml
from utils_io import (
    print_header, print_section,
    safe_input, menu_select,
    write_summary_report, confirm
)
from utils_pluto import (
    load_local_data, load_pluto_data,
    detect_factors, generate_comparisons
)


def validate_factor_combination(metadata_df, factor1, value1, factor2=None, value2=None):
    """Validate that factor combinations exist in the metadata"""
    
    def safe_compare(series, target_value):
        """Safe comparison that handles string/numeric type mismatches"""
        try:
            # Try direct comparison first
            mask1 = series == target_value
            if mask1.sum() > 0:
                return mask1
            
            # If no matches, try converting target_value to match series dtype
            if series.dtype == 'object' or 'str' in str(series.dtype):
                # Series is string, try converting target to string
                mask2 = series == str(target_value)
                if mask2.sum() > 0:
                    return mask2
            
            # If series is numeric, try converting target to numeric
            if 'int' in str(series.dtype) or 'float' in str(series.dtype):
                try:
                    numeric_target = int(target_value) if str(target_value).isdigit() else float(target_value)
                    mask3 = series == numeric_target
                    if mask3.sum() > 0:
                        return mask3
                except (ValueError, TypeError):
                    pass
            
            # Return empty mask if no matches
            return series == target_value  # Original comparison (will be False for all)
            
        except Exception as e:
            print(f"⚠ Comparison error for {target_value}: {e}")
            return series == target_value
    
    if factor2 is None:
        # Single factor validation
        mask = safe_compare(metadata_df[factor1], value1)
        count = mask.sum()
        if count == 0:
            print(f"⚠ Warning: No samples found for {factor1}={value1}")
            print(f"   Available {factor1} values: {sorted(metadata_df[factor1].unique())}")
            return False, 0
        print(f"✓ Found {count} samples for {factor1}={value1}")
        return True, count
    else:
        # Two factor validation
        mask1 = safe_compare(metadata_df[factor1], value1)
        mask2 = safe_compare(metadata_df[factor2], value2)
        mask = mask1 & mask2
        count = mask.sum()
        if count == 0:
            print(f"⚠ Warning: No samples found for {factor1}={value1} AND {factor2}={value2}")
            return False, 0
        print(f"✓ Found {count} samples for {factor1}={value1} AND {factor2}={value2}")
        return True, count


def create_contrast_name_mapping(contrast_name):
    """
    Create a human-readable display name and safe file name from contrast.
    Returns: (original, display, safe_file_name)
    
    Example:
        Input:  "donor:19503_shear_stress:High (15 dynes/cm2)_vs_Static (0 dynes/cm2)"
        Output: 
            original: "donor:19503_shear_stress:High (15 dynes/cm2)_vs_Static (0 dynes/cm2)"
            display:  "donor 19503 shear_stress High (15 dynes/cm2) vs Static (0 dynes/cm2)"
            safe:     "donor_19503_shear_stress_High_15_dynes_cm2_vs_Static_0_dynes_cm2"
    """
    # Display name (human-readable)
    display = contrast_name.replace(":", " ").replace("_vs_", " vs ")
    
    # Safe file name (for directories, no special chars)
    safe = contrast_name
    safe = re.sub(r'\([^)]*\)', lambda m: m.group(0).replace('(', '').replace(')', '').replace('/', '_').replace(' ', '_'), safe)
    safe = re.sub(r'[/:,\s()]', '_', safe)
    safe = re.sub(r'_+', '_', safe)
    safe = re.sub(r'^_|_$', '', safe)
    
    return contrast_name, display, safe


def define_universal_custom_contrasts(metadata_df, factor_info):
    """Universal custom contrast definition that works with any dataset"""
    from utils_io import menu_select, safe_input
    
    print(f"\n📋 Custom Contrast Definition")
    print("=" * 60)
    
    contrast_types = [
        "Between Subject (compare different subjects/donors)",
        "Within Subject (same subject, different conditions)",
        "Pairwise Comparison (any two groups)"
    ]
    
    all_contrasts = []
    
    while True:
        contrast_type = menu_select(contrast_types, "Select contrast type")
        
        if contrast_type == "Between Subject (compare different subjects/donors)":
            new_contrasts = define_single_between_subject_contrast(metadata_df, factor_info)
        elif contrast_type == "Within Subject (same subject, different conditions)":
            new_contrasts = define_single_within_subject_contrast(metadata_df, factor_info)
        else:  # Pairwise Comparison
            new_contrasts = define_single_pairwise_contrast(metadata_df, factor_info)
        
        if new_contrasts:
            all_contrasts.extend(new_contrasts)
        
        # Show current contrasts with formatting
        if all_contrasts:
            print("\n" + "─" * 60)
            print("Current Contrasts:")
            print("─" * 60)
            for i, contrast in enumerate(all_contrasts, 1):
                original, display, safe = create_contrast_name_mapping(contrast)
                print(f"  {i}. {display}")
                print(f"     → File: {safe}/")
            print("─" * 60)
            
            continue_options = [
                "Yes, create another contrast in same type",
                "No, continue with current contrasts",
                "Select another contrast type"
            ]
            
            next_action = menu_select(continue_options, "Create another contrast?")
            
            if next_action == "No, continue with current contrasts":
                break
            elif next_action == "Yes, create another contrast in same type":
                continue
            else:  # "Select another contrast type"
                continue
    
    return all_contrasts


def define_single_between_subject_contrast(metadata_df, factor_info):
    """Define a single between-subject contrast"""
    from utils_io import menu_select, safe_input
    
    print(f"\n🔬 Between Subject Contrasts")
    print("Example: Donor 19503 vs Donor 21452 (with High shear stress)")
    print(f"Available factors: {list(factor_info.keys())}")
    
    # Select subject factor
    subject_factors = [f for f in factor_info.keys() 
                      if 'donor' in f.lower() or 'patient' in f.lower() or 'subject' in f.lower()]
    
    if not subject_factors:
        subject_factors = list(factor_info.keys())
    
    if len(subject_factors) == 1:
        subject_factor = subject_factors[0]
        print(f"Using subject factor: {subject_factor}")
    else:
        subject_factor = menu_select(subject_factors, "Select subject factor")
    
    subject_levels = factor_info[subject_factor]
    if len(subject_levels) < 2:
        print(f"⚠ Need at least 2 {subject_factor} levels for comparison")
        return []
    
    # Select two subjects to compare
    subject1 = menu_select(subject_levels, "Select first subject")
    remaining_subjects = [s for s in subject_levels if s != subject1]
    subject2 = menu_select(remaining_subjects, "Select second subject")
    
    # Validate sample existence
    valid1, count1 = validate_factor_combination(metadata_df, subject_factor, subject1)
    valid2, count2 = validate_factor_combination(metadata_df, subject_factor, subject2)
    
    if not (valid1 and valid2):
        print("⚠ Cannot create contrast - missing samples")
        return []
    
    # Optionally specify condition context
    other_factors = [f for f in factor_info.keys() if f != subject_factor]
    
    if other_factors:
        add_context = menu_select(["Yes", "No"], "Add condition context (e.g., 'High shear stress')?")
        
        if add_context == "Yes":
            context_factor = menu_select(other_factors, "Select context factor")
            context_level = menu_select(factor_info[context_factor], f"Select {context_factor} level")
            
            # Validate combination
            valid1_ctx, _ = validate_factor_combination(metadata_df, subject_factor, subject1, context_factor, context_level)
            valid2_ctx, _ = validate_factor_combination(metadata_df, subject_factor, subject2, context_factor, context_level)
            
            if not (valid1_ctx and valid2_ctx):
                print("⚠ Cannot create contrast - factor combination doesn't exist in data")
                return []
            
            contrast_name = f"{subject_factor}:{subject1}_vs_{subject2}_{context_factor}:{context_level}"
            print(f"✓ Generated contrast: {subject1} vs {subject2} (in {context_level})")
        else:
            contrast_name = f"{subject_factor}:{subject1}_vs_{subject2}"
            print(f"✓ Generated contrast: {subject1} vs {subject2}")
    else:
        contrast_name = f"{subject_factor}:{subject1}_vs_{subject2}"
        print(f"✓ Generated contrast: {subject1} vs {subject2}")
    
    return [contrast_name]


def define_single_within_subject_contrast(metadata_df, factor_info):
    """Define a single within-subject contrast"""
    from utils_io import menu_select, safe_input
    
    print(f"\n🔄 Within Subject Contrasts")
    print("Example: 19503 High vs Static (same donor, different conditions)")
    print(f"Available factors: {list(factor_info.keys())}")
    
    available_factors = list(factor_info.keys())
    if len(available_factors) < 2:
        print("⚠ Need at least 2 factors for within-subject comparison")
        return []
    
    subject_factor = menu_select(available_factors, "Select subject factor (e.g., donor, patient)")
    
    # Select specific subject
    subject_levels = factor_info[subject_factor]
    selected_subject = menu_select(subject_levels, f"Select {subject_factor}")
    
    # Validate subject exists
    valid_subject, _ = validate_factor_combination(metadata_df, subject_factor, selected_subject)
    if not valid_subject:
        return []
    
    # Select condition factor
    condition_factors = [f for f in factor_info.keys() if f != subject_factor]
    if not condition_factors:
        print("⚠ No condition factors available for within-subject comparison")
        return []
    
    condition_factor = menu_select(condition_factors, "Select condition factor")
    condition_levels = factor_info[condition_factor]
    
    if len(condition_levels) < 2:
        print(f"⚠ Need at least 2 {condition_factor} levels for comparison")
        return []
    
    # Select two conditions
    condition1 = menu_select(condition_levels, "Select first condition")
    remaining_conditions = [c for c in condition_levels if c != condition1]
    condition2 = menu_select(remaining_conditions, "Select second condition")
    
    # Validate combinations
    valid1, _ = validate_factor_combination(metadata_df, subject_factor, selected_subject, condition_factor, condition1)
    valid2, _ = validate_factor_combination(metadata_df, subject_factor, selected_subject, condition_factor, condition2)
    
    if not (valid1 and valid2):
        print("⚠ Cannot create contrast - factor combination doesn't exist in data")
        return []
    
    contrast_name = f"{subject_factor}:{selected_subject}_{condition_factor}:{condition1}_vs_{condition2}"
    print(f"✓ Generated contrast: {selected_subject} - {condition1} vs {condition2}")
    
    return [contrast_name]


def define_single_pairwise_contrast(metadata_df, factor_info):
    """Define a single pairwise contrast"""
    from utils_io import menu_select, safe_input
    
    print(f"\n⚖️ Pairwise Contrasts")
    print("Example: Plate1 vs Plate4")
    print(f"Available factors: {list(factor_info.keys())}")
    
    available_factors = list(factor_info.keys())
    selected_factor = menu_select(available_factors, "Select factor for comparison")
    
    factor_levels = factor_info[selected_factor]
    if len(factor_levels) < 2:
        print(f"⚠ Need at least 2 {selected_factor} levels for comparison")
        return []
    
    # Select two levels to compare
    level1 = menu_select(factor_levels, "Select first level")
    remaining_levels = [l for l in factor_levels if l != level1]
    level2 = menu_select(remaining_levels, "Select second level")
    
    # Validate
    valid1, _ = validate_factor_combination(metadata_df, selected_factor, level1)
    valid2, _ = validate_factor_combination(metadata_df, selected_factor, level2)
    
    if not (valid1 and valid2):
        return []
    
    contrast_name = f"{selected_factor}:{level1}_vs_{level2}"
    print(f"✓ Generated contrast: {level1} vs {level2}")
    
    return [contrast_name]


def create_results_structure(project_name, selected_analyses, comparisons):
    """Create results folder structure: Results/ExperimentID/AnalysisName/ContrastName/"""
    import os
    
    base_dir = f"results/{project_name}"
    os.makedirs(base_dir, exist_ok=True)
    
    for analysis in selected_analyses:
        analysis_dir = f"{base_dir}/{analysis.replace(' ', '_')}"
        os.makedirs(analysis_dir, exist_ok=True)
        
        # Create contrast-specific folders for DE
        if analysis == "Differential Expression" and comparisons:
            for contrast in comparisons:
                # Use safe file name for directory
                _, _, safe_name = create_contrast_name_mapping(contrast)
                contrast_dir = f"{analysis_dir}/{safe_name}"
                os.makedirs(contrast_dir, exist_ok=True)
    
    print(f"✓ Results structure created in: {base_dir}/")
    return base_dir


def main():
    # Import needed functions
    from utils_io import safe_input
    from utils_pluto import get_api_token
    
    print_header("RNA-seq Downstream Analysis Setup Wizard")
    print("ℹ️  Tip: Type 'exit' at any prompt to quit the setup wizard")
    print("═" * 60)
    
    # --------------------------------------------------------
    # Step 1 — Data Source Selection
    # --------------------------------------------------------
    print_section("Choose data source")
    mode = menu_select(
        ["Pluto Experiment ID", "Local CSV Files"],
        "Select an input mode"
    )
    
    if mode == "Pluto Experiment ID":
        experiment_id = safe_input("Enter Pluto Experiment ID (e.g. PLX229): ")
        api_token = get_api_token()
        counts_df, metadata_df, counts_path, metadata_path = load_pluto_data(experiment_id, api_token)
        project_name = experiment_id
    else:
        project_name = safe_input("Enter project name (or 'exit' to quit): ")
        if project_name and project_name.lower() == 'exit':
            print("ℹ Exiting setup...")
            return
        
        print("\nEnter file paths (use absolute paths, or 'exit' to quit):")
        counts_path = safe_input("Path to counts CSV: ")
        if counts_path and counts_path.lower() == 'exit':
            print("ℹ Exiting setup...")
            return
        
        metadata_path = safe_input("Path to metadata CSV: ")
        if metadata_path and metadata_path.lower() == 'exit':
            print("ℹ Exiting setup...")
            return
        
        try:
            counts_df, metadata_df = load_local_data(counts_path, metadata_path)
        except FileNotFoundError as e:
            print(f"✗ File not found: {e}")
            print("Please provide absolute paths to your CSV files.")
            return
        except Exception as e:
            print(f"✗ Error loading files: {e}")
            return
    
    # --------------------------------------------------------
    # Step 2 — Analysis Selection
    # --------------------------------------------------------
    print_section("Select Analysis Types")
    analysis_options = [
        "Batch Correction",
        "Power Analysis",
        "Differential Expression",
        "Pathway Analysis",
        "All Analyses",
        "Exit Setup"
    ]
    
    selected_analyses = []
    print("\nSelect analyses to perform:")
    print("You can enter multiple numbers separated by commas (e.g., 1,2,3)")
    
    for i, option in enumerate(analysis_options, 1):
        print(f" {i}. {option}")
    
    while True:
        choice_input = input("→ Select option(s): ").strip()
        if not choice_input:
            continue
        
        try:
            if ',' in choice_input:
                choices = [int(x.strip()) for x in choice_input.split(',')]
                for choice_num in choices:
                    if 1 <= choice_num <= len(analysis_options):
                        selected_option = analysis_options[choice_num - 1]
                        if selected_option == "Exit Setup":
                            if not selected_analyses:
                                print("⚠ No analyses selected. Exiting...")
                                return
                            break
                        elif selected_option == "All Analyses":
                            selected_analyses = ["Batch Correction", "Power Analysis", 
                                               "Differential Expression", "Pathway Analysis"]
                            break
                        elif selected_option not in selected_analyses:
                            selected_analyses.append(selected_option)
                            print(f"✓ Added: {selected_option}")
                break
            else:
                choice_num = int(choice_input)
                if 1 <= choice_num <= len(analysis_options):
                    selected_option = analysis_options[choice_num - 1]
                    if selected_option == "Exit Setup":
                        if not selected_analyses:
                            print("⚠ No analyses selected. Exiting...")
                            return
                        break
                    elif selected_option == "All Analyses":
                        selected_analyses = ["Batch Correction", "Power Analysis", 
                                           "Differential Expression", "Pathway Analysis"]
                        break
                    elif selected_option not in selected_analyses:
                        selected_analyses.append(selected_option)
                        print(f"✓ Added: {selected_option}")
                    
                    add_more = safe_input("Add another analysis? (y/n)", "n")
                    if add_more.lower() not in ['y', 'yes']:
                        break
                else:
                    print("✗ Invalid choice. Try again.")
        except ValueError:
            print("✗ Invalid input. Enter numbers only (e.g., 1 or 1,2,3).")
    
    print(f"\n✓ Selected analyses: {', '.join(selected_analyses)}")
    
    # Initialize configuration variables
    de_tool = None
    filtering_config = {}
    normalization_config = {}
    design_complexity = "N/A"
    applied_rule = "N/A"
    
    # Direct DE tool selection if differential expression is selected
    if "Differential Expression" in selected_analyses:
        print_section("Differential Expression Configuration")
        de_tool = menu_select(["DESeq2", "limma"], "Choose DE tool")
        print(f"✓ Selected DE tool: {de_tool}")
    
    # --------------------------------------------------------
    # Step 3 — Infer metadata factors
    # --------------------------------------------------------
    print_section("Detecting factors from metadata")
    factor_candidates = detect_factors(metadata_df)
    
    print("Detected categorical fields with levels:")
    factor_info = {}
    
    for col in factor_candidates:
        unique_values = metadata_df[col].unique()
        levels = [str(x) for x in unique_values if str(x) != 'nan']
        
        try:
            levels = sorted(levels, key=lambda x: (x.replace('.','').replace('-','').isdigit() == False, x))
        except:
            levels = sorted(levels)
        
        factor_info[col] = levels
        sample_counts = metadata_df[col].value_counts().to_dict()
        
        print(f" → {col}: {levels}")
        for level in levels:
            original_level = level
            if level.replace('.','').replace('-','').isdigit():
                try:
                    original_level = int(level) if '.' not in level else float(level)
                except:
                    pass
            count = sample_counts.get(original_level, sample_counts.get(level, 0))
            print(f"    - {level}: {count} samples")
    
    primary_factor = menu_select(factor_candidates, "Select primary factor")
    print(f"Primary factor levels: {factor_info[primary_factor]}")
    
    # Calculate smallest group size for filtering configuration
    if "Differential Expression" in selected_analyses:
        actual_group_sizes = metadata_df[primary_factor].value_counts()
        actual_smallest_group = actual_group_sizes.min()
        actual_total_samples = len(metadata_df)
        
        print("Analyzing dataset characteristics...")
        print(f"   Total samples: {actual_total_samples}")
        print(f"   Smallest group size: {actual_smallest_group}")
        print(f"   Potential factors: {len(factor_candidates)}")
        
        # Configure gene filtering and normalization
        print_section("Gene Filtering and Normalisation")
        
        if de_tool == "DESeq2":
            print("type: raw_count filtering")
            print('deseq2_note: "No normalization - automatic"')
            print("threshold: 10")
            print('logic: "rowSums >= 10"')
            
            filtering_config = {
                "method": "DESeq2",
                "type": "raw_count",
                "threshold": 10,
                "logic": "rowSums >= 10",
                "description": "Keep genes with raw count sum >= 10 across all samples"
            }
            normalization_config = {
                "method": "automatic",
                "note": "No normalization - automatic"
            }
        else:  # limma
            print("limma_voom:")
            print("normalization:")
            print("  current: CPM")
            print("  future: [TMM, RLE, upperquartile]")
            print(f"threshold: 10")
            print(f"min_samples: {actual_smallest_group}")
            print(f'logic: "CPM >= 10 in >= {actual_smallest_group} samples"')
            
            filtering_config = {
                "method": "limma-voom",
                "type": "cpm",
                "threshold": 10,
                "min_samples": actual_smallest_group,
                "logic": f"CPM >= 10 in >= {actual_smallest_group} samples",
                "description": f"Keep genes with CPM >= 10 in at least {actual_smallest_group} samples"
            }
            normalization_config = {
                "method": "CPM",
                "note": "Current: CPM normalization. Future: TMM, RLE, upperquartile",
                "future_options": ["TMM", "RLE", "upperquartile"]
            }
        
        # Ask user for confirmation/customization
        use_default = menu_select(
            ["Yes, use suggested thresholds", "No, let me customize"],
            "Use suggested filtering thresholds?"
        )
        
        if use_default == "No, let me customize":
            if de_tool == "DESeq2":
                threshold_input = safe_input("Enter custom threshold for rowSums", "10")
                try:
                    filtering_config["threshold"] = int(threshold_input)
                    filtering_config["logic"] = f"rowSums >= {filtering_config['threshold']}"
                except ValueError:
                    print("Invalid threshold, using default: 10")
            else:
                cpm_input = safe_input("Enter CPM threshold", "10")
                try:
                    filtering_config["threshold"] = float(cpm_input)
                    filtering_config["logic"] = f"CPM >= {filtering_config['threshold']} in >= {actual_smallest_group} samples"
                except ValueError:
                    print("Invalid input, using defaults")
    
    batch_factor = menu_select(["None"] + factor_candidates, "Select batch factor (or None)")
    batch_factor = None if batch_factor == "None" else batch_factor
    
    if batch_factor:
        print(f"Batch factor levels: {factor_info[batch_factor]}")
        
        # Check ComBat-seq batch size requirements
        batch_table = metadata_df[batch_factor].value_counts()
        min_samples_per_batch = batch_table.min()
        single_sample_batches = batch_table[batch_table == 1].index.tolist()
        
        if min_samples_per_batch < 2:
            print("\n⚠️  BATCH CORRECTION LIMITATION:")
            print(f"   ComBat-seq requires ≥2 samples per batch (found: {min_samples_per_batch} minimum)")
            if single_sample_batches:
                print(f"   Batches with single samples: {', '.join(map(str, single_sample_batches))}")
            print("   ┌─ Pipeline Decision:")
            print("   │  ✗ Skip ComBat-seq batch correction")
            print("   │  ✓ Use gene-filtered counts instead")
            print("   └─ Power analysis will use filtered counts")
            
            if "Batch Correction" in selected_analyses:
                proceed = confirm("Continue with pipeline using filtered counts (no batch correction)?")
                if not proceed:
                    print("   Setup cancelled. Consider merging batches or using different batch factor.")
                    exit(1)
            else:
                print("   ✓ Will proceed with gene-filtered counts")
    
    additional_factors = [
        c for c in factor_candidates
        if c not in (primary_factor, batch_factor)
    ]
    
    # --------------------------------------------------------
    # Step 4 — Contrast Definition
    # --------------------------------------------------------
    comparisons = []
    
    if "Differential Expression" in selected_analyses:
        print_section("Define Contrasts for Differential Expression")
        print(f"Primary factor: {primary_factor}")
        comparisons = define_universal_custom_contrasts(metadata_df, factor_info)
        
        # Show final contrast summary
        if comparisons:
            print("\n" + "=" * 60)
            print("FINAL CONTRAST SUMMARY")
            print("=" * 60)
            for i, contrast in enumerate(comparisons, 1):
                original, display, safe = create_contrast_name_mapping(contrast)
                print(f"\n{i}. Display: {display}")
                print(f"   YAML:    {original}")
                print(f"   Folder:  {safe}/")
            print("=" * 60)
    else:
        print("ℹ Skipping contrast definition (Differential Expression not selected)")
    
    # --------------------------------------------------------
    # Step 5 — Build YAML config
    # --------------------------------------------------------
    print_section("Building analysis_config.yaml")
    
    config = {
        "project_name": project_name,
        "input_files": {
            "counts": counts_path,
            "metadata": metadata_path,
        },
        "experimental_design": {
            "primary_factor": primary_factor,
            "batch_factor": batch_factor,
            "additional_factors": additional_factors
        },
        "selected_analyses": selected_analyses,
        "analysis_parameters": {
            "de_tool": de_tool,
            "perform_batch_correction": "Batch Correction" in selected_analyses,
            "perform_power_analysis": "Power Analysis" in selected_analyses,
            "perform_differential_expression": "Differential Expression" in selected_analyses,
            "perform_pathway_analysis": "Pathway Analysis" in selected_analyses,
            "metadata_columns": {
                "condition": primary_factor,
                "batch": batch_factor if batch_factor else None
            },
            "results_structure": f"data/{project_name}/analysis_results"
        },
        "filtering_config": filtering_config,
        "normalization_config": normalization_config,
        "dataset_characteristics": {
            "total_samples": len(metadata_df),
            "smallest_group_size": filtering_config.get("min_samples", "N/A"),
            "design_complexity": design_complexity if "Differential Expression" in selected_analyses else "N/A",
            "applied_rule": applied_rule if "Differential Expression" in selected_analyses else "N/A"
        },
        "comparisons": comparisons  # Store original format
    }
    
    with open("analysis_config.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False, default_flow_style=False)
    
    print("✓ analysis_config.yaml saved successfully")
    print("\n" + "=" * 60)
    print("✅ Setup Complete!")
    print("=" * 60)
    print(f"Configuration saved: analysis_config.yaml")
    print(f"Next step: Run Nextflow pipeline")
    print("=" * 60)


if __name__ == "__main__":
    main()