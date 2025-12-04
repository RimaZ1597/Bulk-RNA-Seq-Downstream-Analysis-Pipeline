#!/usr/bin/env python3
"""
Multi-Dataset Interactive Setup for RNA-seq Analysis Pipeline
Supports multiple studies with dataset selection
Enhanced with multi-factor comparison options and selective comparison choice
"""
import os
import json
import yaml
import pandas as pd
from pathlib import Path
from collections import defaultdict
from itertools import product

def detect_multi_dataset_structure():
    """Detect multiple dataset structure"""
    print("=== DETECTING MULTI-DATASET STRUCTURE ===")
    data_info = {
        'studies': {},
        'datasets_found': 0
    }
    data_dir = Path("data")
    if not data_dir.exists():
        print("❌ Data directory not found")
        return data_info
    
    # Look for study directories
    study_dirs = [d for d in data_dir.iterdir() if d.is_dir()]
    for study_dir in study_dirs:
        study_name = study_dir.name
        print(f"\n📁 Checking study: {study_name}")
        study_info = {
            'path': str(study_dir),
            'count_files': [],
            'metadata_files': [],
            'cleaned_metadata': None,
            'valid': False
        }
        
        # Look for count files
        count_patterns = ['*.csv', '*.tsv', '*.txt']
        for pattern in count_patterns:
            count_files = list(study_dir.glob(f"*count*{pattern[1:]}")) + \
                         list(study_dir.glob(f"count*{pattern[1:]}")) + \
                         list(study_dir.glob(f"*expression*{pattern[1:]}"))
            study_info['count_files'].extend([str(f) for f in count_files])
        
        # Look for metadata files  
        for pattern in count_patterns:
            meta_files = list(study_dir.glob(f"*meta*{pattern[1:]}")) + \
                        list(study_dir.glob(f"*sample*{pattern[1:]}")) + \
                        list(study_dir.glob(f"metadata{pattern[1:]}"))
            study_info['metadata_files'].extend([str(f) for f in meta_files])
        
        # Look for cleaned metadata
        cleaned_dir = study_dir / "cleaned_metadata"
        if cleaned_dir.exists():
            cleaned_files = list(cleaned_dir.glob("*.csv"))
            if cleaned_files:
                study_info['cleaned_metadata'] = str(cleaned_files[0])
        
        # Check if study is valid (has both count and metadata)
        if study_info['count_files'] and study_info['metadata_files']:
            study_info['valid'] = True
            data_info['datasets_found'] += 1
            print(f"   ✅ Valid study found:")
            print(f"      Count files: {len(study_info['count_files'])}")
            print(f"      Metadata files: {len(study_info['metadata_files'])}")
            if study_info['cleaned_metadata']:
                print(f"      Cleaned metadata: ✅")
        else:
            print(f"   ❌ Invalid study (missing count or metadata files)")
        
        data_info['studies'][study_name] = study_info
    
    return data_info

def select_dataset(data_info):
    """Let user select which dataset to analyze"""
    valid_studies = {name: info for name, info in data_info['studies'].items() if info['valid']}
    if not valid_studies:
        print("❌ No valid datasets found")
        return None
    
    print(f"\n=== DATASET SELECTION ===")
    print(f"Found {len(valid_studies)} valid dataset(s):")
    
    study_list = list(valid_studies.keys())
    for i, study_name in enumerate(study_list):
        info = valid_studies[study_name]
        print(f"   {i+1}) {study_name}")
        print(f"      📁 Path: {info['path']}")
        print(f"      📊 Count files: {len(info['count_files'])}")
        print(f"      📋 Metadata files: {len(info['metadata_files'])}")
        if info['cleaned_metadata']:
            print(f"      🧹 Cleaned metadata: Available")
        print()
    
    if len(study_list) == 1:
        selected_study = study_list[0]
        print(f"✅ Using dataset: {selected_study}")
    else:
        while True:
            try:
                choice = input(f"Select dataset (1-{len(study_list)}): ")
                choice_idx = int(choice) - 1
                if 0 <= choice_idx < len(study_list):
                    selected_study = study_list[choice_idx]
                    print(f"✅ Selected dataset: {selected_study}")
                    break
                else:
                    print("❌ Invalid choice, please try again")
            except ValueError:
                print("❌ Please enter a number")
    
    return selected_study, valid_studies[selected_study]

def select_files_for_study(study_info):
    """Select specific files within a study"""
    selected_files = {}
    
    # Select count file
    count_files = study_info['count_files']
    if len(count_files) == 1:
        selected_files['counts'] = count_files[0]
        print(f"✅ Using count file: {Path(selected_files['counts']).name}")
    else:
        print("\nAvailable count files:")
        for i, f in enumerate(count_files):
            print(f"   {i+1}) {Path(f).name}")
        choice = int(input("Select count file (number): ")) - 1
        selected_files['counts'] = count_files[choice]
    
    # Select metadata file (prioritize cleaned if available)
    if study_info['cleaned_metadata']:
        print(f"\nFound cleaned metadata: {Path(study_info['cleaned_metadata']).name}")
        use_cleaned = input("Use cleaned metadata? (y/n): ")
        if use_cleaned.lower() == 'y':
            selected_files['metadata'] = study_info['cleaned_metadata']
            selected_files['metadata_type'] = 'cleaned'
            print("✅ Using cleaned metadata")
        else:
            selected_files['metadata'] = select_original_metadata(study_info['metadata_files'])
            selected_files['metadata_type'] = 'original'
    else:
        selected_files['metadata'] = select_original_metadata(study_info['metadata_files'])
        selected_files['metadata_type'] = 'original'
    
    return selected_files

def select_original_metadata(metadata_files):
    """Select from original metadata files"""
    if len(metadata_files) == 1:
        selected_file = metadata_files[0]
        print(f"✅ Using metadata file: {Path(selected_file).name}")
    else:
        print("\nAvailable metadata files:")
        for i, f in enumerate(metadata_files):
            print(f"   {i+1}) {Path(f).name}")
        choice = int(input("Select metadata file (number): ")) - 1
        selected_file = metadata_files[choice]
    
    return selected_file

def validate_sample_matching(count_file, metadata_file, sample_id_column):
    """Validate that metadata samples match count matrix columns"""
    print("\n=== VALIDATING SAMPLE MATCHING ===")
    
    # Read count matrix headers
    count_df = pd.read_csv(count_file, nrows=0)  # Just headers
    count_samples = list(count_df.columns[1:])  # Skip first column (gene names)
    
    # Read metadata samples
    meta_df = pd.read_csv(metadata_file)
    meta_samples = list(meta_df[sample_id_column])
    
    # Check overlap
    common_samples = set(count_samples) & set(meta_samples)
    missing_in_counts = set(meta_samples) - set(count_samples)
    missing_in_metadata = set(count_samples) - set(meta_samples)
    
    print(f"✅ Count matrix samples: {len(count_samples)}")
    print(f"✅ Metadata samples: {len(meta_samples)}")
    print(f"✅ Common samples: {len(common_samples)}")
    
    if missing_in_counts:
        print(f"⚠️  Samples in metadata but not in counts: {len(missing_in_counts)}")
    if missing_in_metadata:
        print(f"⚠️  Samples in counts but not in metadata: {len(missing_in_metadata)}")
    
    if len(common_samples) == 0:
        print("❌ No matching samples found!")
        return False
    
    return True

def analyze_metadata(metadata_file):
    """Analyze metadata structure and detect experimental factors"""
    print(f"\n=== ANALYZING METADATA STRUCTURE ===")
    
    # Try different separators
    separators = [',', '\t', ';']
    df = None
    for sep in separators:
        try:
            df = pd.read_csv(metadata_file, sep=sep)
            if len(df.columns) > 1:
                break
        except:
            continue
    
    if df is None:
        print("❌ Could not read metadata file")
        return None, None
    
    print(f"✅ Detected {len(df.columns)} columns in metadata")
    print(f"✅ Loaded {len(df)} samples")
    
    # Show column information
    metadata_info = {}
    print("\nColumn analysis:")
    
    for col in df.columns:
        unique_vals = df[col].nunique()
        sample_vals = df[col].dropna().unique()[:5]  # Show first 5 unique values
        
        print(f"   {col}: {unique_vals} unique values")
        print(f"      Sample values: {list(sample_vals)}")
        
        # Categorize column type
        if unique_vals == len(df):
            col_type = 'identifier'
        elif unique_vals <= 10:
            col_type = 'categorical'
        elif df[col].dtype in ['int64', 'float64']:
            col_type = 'numeric'
        else:
            col_type = 'other'
        
        metadata_info[col] = {
            'type': col_type,
            'unique_count': unique_vals,
            'sample_values': [str(val) if hasattr(val, 'dtype') else val for val in sample_vals]
        }
    
    return metadata_info, df

def interactive_column_deletion(metadata_df, metadata_info):
    """Interactive column deletion functionality"""
    print("\n=== COLUMN DELETION OPTIONS ===")
    
    current_columns = list(metadata_df.columns)
    columns_to_delete = []
    
    print(f"\nCurrent metadata has {len(current_columns)} columns:")
    for i, col in enumerate(current_columns):
        col_info = metadata_info[col]
        print(f"   {i+1:2d}) {col:20s} ({col_info['type']:12s}) - {col_info['unique_count']} unique values")
    
    print("\nColumn deletion options:")
    print("   1) Select specific columns to delete")
    print("   2) Keep only selected columns (delete the rest)")
    print("   3) Delete by column type")
    print("   4) No deletion - keep all columns")
    print("   0) Exit setup")
    
    choice = input("\nSelect deletion method (1/2/3/4/0): ")
    
    if choice == '0':
        print("🚪 Exiting setup...")
        exit(0)
    
    if choice == '1':
        # Select specific columns to delete
        print("\nSelect columns to DELETE:")
        print("Enter column numbers separated by commas (e.g., 1,3,5)")
        print("Or enter column names separated by commas")
        print("Or enter 'none' to skip deletion")
        
        delete_input = input("Columns to delete: ").strip()
        
        if delete_input.lower() != 'none':
            if delete_input.replace(',', '').replace(' ', '').isdigit():
                # Numbers provided
                indices = [int(x.strip()) - 1 for x in delete_input.split(',')]
                columns_to_delete = [current_columns[i] for i in indices if 0 <= i < len(current_columns)]
            else:
                # Column names provided
                columns_to_delete = [col.strip() for col in delete_input.split(',') if col.strip() in current_columns]
    
    elif choice == '2':
        # Keep only selected columns
        print("\nSelect columns to KEEP:")
        print("Enter column numbers separated by commas (e.g., 1,2,4)")
        print("Or enter column names separated by commas")
        
        keep_input = input("Columns to keep: ").strip()
        
        if keep_input.replace(',', '').replace(' ', '').isdigit():
            # Numbers provided
            indices = [int(x.strip()) - 1 for x in keep_input.split(',')]
            columns_to_keep = [current_columns[i] for i in indices if 0 <= i < len(current_columns)]
        else:
            # Column names provided
            columns_to_keep = [col.strip() for col in keep_input.split(',') if col.strip() in current_columns]
        
        columns_to_delete = [col for col in current_columns if col not in columns_to_keep]
    
    elif choice == '3':
        # Delete by column type
        print("\nDelete columns by type:")
        print("   1) Delete all identifier columns")
        print("   2) Delete all categorical columns")
        print("   3) Delete all numeric columns")
        print("   4) Delete all 'other' type columns")
        print("   5) Custom type selection")
        
        type_choice = input("Select type to delete (1/2/3/4/5): ")
        
        type_mapping = {
            '1': 'identifier',
            '2': 'categorical',
            '3': 'numeric',
            '4': 'other'
        }
        
        if type_choice in type_mapping:
            target_type = type_mapping[type_choice]
            columns_to_delete = [col for col, info in metadata_info.items() if info['type'] == target_type]
        elif type_choice == '5':
            print("Available types:", set(info['type'] for info in metadata_info.values()))
            types_to_delete = input("Enter types to delete (comma-separated): ").split(',')
            types_to_delete = [t.strip() for t in types_to_delete]
            columns_to_delete = [col for col, info in metadata_info.items() if info['type'] in types_to_delete]
    
    elif choice == '4':
        # No deletion
        print("✅ Keeping all columns")
        columns_to_delete = []
    
    # Confirm deletion
    if columns_to_delete:
        print(f"\n📋 Columns selected for deletion ({len(columns_to_delete)}):")
        for col in columns_to_delete:
            print(f"   - {col}")
        
        confirm = input(f"\nConfirm deletion of {len(columns_to_delete)} columns? (y/n): ")
        
        if confirm.lower() == 'y':
            # Perform deletion
            cleaned_df = metadata_df.drop(columns=columns_to_delete)
            
            # Update metadata_info
            cleaned_metadata_info = {col: info for col, info in metadata_info.items()
                                   if col not in columns_to_delete}
            
            print(f"✅ Deleted {len(columns_to_delete)} columns")
            print(f"✅ Remaining columns: {len(cleaned_df.columns)}")
            
            # Show remaining columns
            print("\nRemaining columns:")
            for col in cleaned_df.columns:
                print(f"   - {col}")
            
            return cleaned_metadata_info, cleaned_df, columns_to_delete
        else:
            print("❌ Column deletion cancelled")
            return metadata_info, metadata_df, []
    else:
        return metadata_info, metadata_df, []

def save_cleaned_metadata(cleaned_df, original_file, deleted_columns, study_path):
    """Save the cleaned metadata file within study directory"""
    # Create cleaned metadata directory within study
    study_dir = Path(study_path)
    cleaned_dir = study_dir / "cleaned_metadata"
    cleaned_dir.mkdir(exist_ok=True)
    
    # Generate new filename
    original_path = Path(original_file)
    cleaned_filename = cleaned_dir / f"cleaned_{original_path.name}"
    
    # Save cleaned metadata
    cleaned_df.to_csv(cleaned_filename, index=False)
    
    # Save deletion log
    deletion_log = {
        'original_file': original_file,
        'cleaned_file': str(cleaned_filename),
        'deleted_columns': deleted_columns,
        'original_column_count': len(cleaned_df.columns) + len(deleted_columns),
        'final_column_count': len(cleaned_df.columns)
    }
    
    with open(cleaned_dir / "deletion_log.yaml", 'w') as f:
        yaml.safe_dump(deletion_log, f, default_flow_style=False)
    
    print(f"✅ Cleaned metadata saved to: {cleaned_filename}")
    print(f"✅ Deletion log saved to: {cleaned_dir / 'deletion_log.yaml'}")
    
    return str(cleaned_filename)

def get_experimental_design(metadata_info, metadata_df):
    """Interactive selection of experimental design"""
    config = {}
    
    print("\n=== EXPERIMENTAL DESIGN SETUP ===")
    
    # 1. Sample ID column
    print("\n1. SAMPLE IDENTIFIER:")
    id_candidates = [col for col, info in metadata_info.items()
                    if info['type'] == 'identifier' or 'sample' in col.lower() or 'id' in col.lower()]
    
    if id_candidates:
        print("Suggested sample ID columns:")
        for i, col in enumerate(id_candidates):
            print(f"   {i+1}) {col}")
        print(f"   0) Other column")
        
        choice = input("Select sample ID column (number or column name): ")
        if choice.isdigit() and int(choice) > 0 and int(choice) <= len(id_candidates):
            config['sample_id_column'] = id_candidates[int(choice)-1]
        elif choice == '0':
            available_cols = list(metadata_info.keys())
            print("Available columns:", available_cols)
            config['sample_id_column'] = input("Enter sample ID column name: ")
        else:
            config['sample_id_column'] = choice
    else:
        print("Available columns:")
        for col in metadata_info.keys():
            print(f"   {col}")
        config['sample_id_column'] = input("Enter sample ID column name: ")
    
    # 2. Treatment/Condition factors
    print("\n2. EXPERIMENTAL FACTORS:")
    categorical_cols = [col for col, info in metadata_info.items()
                       if info['type'] == 'categorical' and col != config['sample_id_column']]
    
    print("Available categorical columns for experimental factors:")
    for i, col in enumerate(categorical_cols):
        values = metadata_info[col]['sample_values']
        print(f"   {i+1}) {col}: {values}")
    
    # Primary factor (main condition)
    print("\n   Select PRIMARY experimental factor (main condition to compare):")
    choice = input("Enter column name or number: ")
    if choice.isdigit() and int(choice) <= len(categorical_cols):
        config['primary_factor'] = categorical_cols[int(choice)-1]
    else:
        config['primary_factor'] = choice
    
    # Additional factors
    config['additional_factors'] = []
    remaining_cols = [col for col in categorical_cols if col != config['primary_factor']]
    
    if remaining_cols:
        print(f"\n   Available additional factors: {remaining_cols}")
        print("   Select ADDITIONAL factors (batch, treatment, time, etc.):")
        print("   Enter column names separated by commas, or 'none':")
        
        additional = input("Additional factors: ").strip()
        if additional.lower() != 'none':
            if additional.replace(',', '').replace(' ', '').isdigit():
                # Numbers provided
                indices = [int(x.strip()) - 1 for x in additional.split(',')]
                config['additional_factors'] = [remaining_cols[i] for i in indices if 0 <= i < len(remaining_cols)]
            else:
                # Column names provided
                config['additional_factors'] = [f.strip() for f in additional.split(',')]
    
    # 3. Batch correction
    print("\n3. BATCH CORRECTION:")
    all_factors = [config['primary_factor']] + config['additional_factors']
    batch_candidates = [col for col in all_factors
                       if any(term in col.lower() for term in ['batch', 'plate', 'run', 'lane'])]
    
    if batch_candidates:
        print(f"Potential batch factors detected: {batch_candidates}")
        config['batch_factor'] = input("Enter batch factor column name (or 'none'): ")
        if config['batch_factor'].lower() == 'none':
            config['batch_factor'] = None
    else:
        if all_factors:
            print(f"Available factors for batch correction: {all_factors}")
        config['batch_factor'] = input("Enter batch factor column name (or 'none' if no batch effect): ")
        if config['batch_factor'].lower() == 'none':
            config['batch_factor'] = None
    
    return config

def get_analysis_parameters():
    """Get analysis thresholds - simplified version with defaults"""
    print("\n=== ANALYSIS PARAMETERS ===")
    
    params = {}
    
    # Statistical thresholds
    print("\n1. STATISTICAL THRESHOLDS:")
    params['pvalue_threshold'] = float(input("P-value threshold (default 0.05): ") or "0.05")
    params['padj_threshold'] = float(input("Adjusted p-value threshold (default 0.05): ") or "0.05")
    params['logfc_threshold'] = float(input("Log fold-change threshold (default 1.0): ") or "1.0")
    
    # Filtering parameters
    print("\n2. FILTERING PARAMETERS:")
    params['min_counts'] = int(input("Minimum count threshold (default 10): ") or "10")
    
    # Set default analysis options (all enabled, CSV only)
    params['min_samples'] = None
    params['perform_pca'] = False  # Disabled for CSV only
    params['generate_heatmaps'] = False  # Disabled for CSV only
    params['save_normalized_counts'] = True
    params['generate_reports'] = True
    
    print("✅ Analysis options set: CSV only")
    return params

def generate_combination_options(metadata_info, all_factors):
    """Generate all possible factor combination options for selection"""
    options = []
    
    # Get all factor combinations
    factor_combinations = {}
    for factor in all_factors:
        factor_combinations[factor] = metadata_info[factor]['sample_values']
    
    # Generate all possible group combinations
    all_combinations = list(product(*[factor_combinations[f] for f in all_factors]))
    
    for combo in all_combinations:
        combo_dict = {all_factors[i]: combo[i] for i in range(len(all_factors))}
        combo_str = "_".join([f"{factor}:{combo_dict[factor]}" for factor in all_factors])
        options.append({
            'string': combo_str,
            'dict': combo_dict,
            'tuple': combo
        })
    
    return options

def select_specific_comparisons(comparisons):
    """Allow user to select specific comparisons from the generated list"""
    print(f"\n=== COMPARISON SELECTION ===")
    print(f"Generated {len(comparisons)} total comparisons:")
    
    # Show all comparisons with numbers
    for i, comp in enumerate(comparisons, 1):
        print(f"   {i:2d}) {comp}")
    
    if len(comparisons) <= 10:
        print("\nSmall number of comparisons - proceeding with all.")
        return comparisons
    
    print(f"\n⚠️  Large number of comparisons ({len(comparisons)}) detected!")
    print("Options:")
    print("   1) Use all comparisons")
    print("   2) Select specific comparisons by number")
    print("   3) Select range of comparisons")
    print("   0) Exit setup")
    
    choice = input("Choose option (1/2/3/0): ").strip()
    
    if choice == '0':
        print("🚪 Exiting setup...")
        exit(0)
    elif choice == '1':
        print("✅ Using all comparisons")
        return comparisons
    elif choice == '2':
        return select_by_numbers(comparisons)
    elif choice == '3':
        return select_by_range(comparisons)
    else:
        print("❌ Invalid choice, using all comparisons")
        return comparisons

def select_by_numbers(comparisons):
    """Select specific comparisons by number"""
    print(f"\n=== SELECT BY NUMBERS ===")
    print("Enter comparison numbers separated by commas")
    print("Examples: 1,8,9,4  or  1,3,5,10,15")
    print("Or enter 'all' to use all comparisons")
    
    while True:
        user_input = input("Enter numbers: ").strip()
        
        if user_input.lower() == 'all':
            return comparisons
        
        try:
            # Parse numbers
            numbers = [int(x.strip()) for x in user_input.split(',')]
            
            # Validate numbers
            invalid_numbers = [n for n in numbers if n < 1 or n > len(comparisons)]
            if invalid_numbers:
                print(f"❌ Invalid numbers: {invalid_numbers}")
                print(f"   Valid range: 1-{len(comparisons)}")
                continue
            
            # Remove duplicates and sort
            numbers = sorted(list(set(numbers)))
            
            # Get selected comparisons
            selected_comparisons = [comparisons[n-1] for n in numbers]
            
            print(f"\n📋 Selected {len(selected_comparisons)} comparisons:")
            for i, comp in enumerate(selected_comparisons, 1):
                print(f"   {i}) {comp}")
            
            confirm = input(f"\nConfirm selection of {len(selected_comparisons)} comparisons? (y/n): ")
            if confirm.lower() == 'y':
                return selected_comparisons
            else:
                print("Let's try again...")
                continue
                
        except ValueError:
            print("❌ Invalid input. Please enter numbers separated by commas (e.g., 1,8,9,4)")
            continue

def select_by_range(comparisons):
    """Select comparisons by range"""
    print(f"\n=== SELECT BY RANGE ===")
    print("Enter start and end numbers")
    print("Examples: 1-10, 5-15, 20-30")
    
    while True:
        user_input = input("Enter range (start-end): ").strip()
        
        try:
            if '-' not in user_input:
                print("❌ Please use format: start-end (e.g., 1-10)")
                continue
            
            start, end = user_input.split('-')
            start = int(start.strip())
            end = int(end.strip())
            
            # Validate range
            if start < 1 or end > len(comparisons) or start > end:
                print(f"❌ Invalid range. Valid range: 1-{len(comparisons)}")
                continue
            
            # Get selected comparisons
            selected_comparisons = comparisons[start-1:end]
            
            print(f"\n📋 Selected comparisons {start}-{end} ({len(selected_comparisons)} total):")
            for i, comp in enumerate(selected_comparisons, start):
                print(f"   {i}) {comp}")
            
            confirm = input(f"\nConfirm selection of {len(selected_comparisons)} comparisons? (y/n): ")
            if confirm.lower() == 'y':
                return selected_comparisons
            else:
                print("Let's try again...")
                continue
                
        except ValueError:
            print("❌ Invalid input. Please use format: start-end (e.g., 1-10)")
            continue

def define_comparisons(metadata_info, exp_design):
    """Define which comparisons to perform - Enhanced for multi-factor comparisons with selection"""
    print("\n=== COMPARISON SETUP ===")
    
    primary_factor = exp_design['primary_factor']
    additional_factors = exp_design['additional_factors']
    all_factors = [primary_factor] + additional_factors
    
    print(f"\nAvailable factors for comparisons:")
    for i, factor in enumerate(all_factors):
        levels = metadata_info[factor]['sample_values']
        print(f"   {i+1}) {factor}: {levels}")
    
    comparisons = []
    
    print("\nComparison options:")
    print("   1) Single factor - Compare levels within primary factor only")
    print("   2) Multi-factor - Compare combinations of all factors")
    print("   3) Custom factor selection - Choose specific factors to combine")
    print("   4) All pairwise comparisons")
    print("   0) Exit setup")
    
    choice = input("Select comparison type (1/2/3/4/0): ")
    
    if choice == '0':
        print("🚪 Exiting setup...")
        exit(0)
    
    if choice == '1':
        # Single factor comparisons
        primary_levels = metadata_info[primary_factor]['sample_values']
        print(f"\nPrimary factor '{primary_factor}' has levels: {primary_levels}")
        
        sub_choice = input("All pairwise (p) or vs control (c)? ")
        if sub_choice.lower() == 'p':
            for i in range(len(primary_levels)):
                for j in range(i+1, len(primary_levels)):
                    comparisons.append(f"{primary_factor}:{primary_levels[j]} vs {primary_factor}:{primary_levels[i]}")
        else:
            control = input(f"Enter control level from {primary_levels}: ")
            for level in primary_levels:
                if level != control:
                    comparisons.append(f"{primary_factor}:{level} vs {primary_factor}:{control}")
    
    elif choice == '2':
        # Multi-factor comparisons
        print(f"\nGenerating all possible combinations from factors: {all_factors}")
        combo_options = generate_combination_options(metadata_info, all_factors)
        
        print(f"\nFound {len(combo_options)} possible group combinations:")
        for i, option in enumerate(combo_options[:10]):  # Show first 10
            print(f"   {i+1:2d}) {option['string']}")
        if len(combo_options) > 10:
            print(f"   ... and {len(combo_options) - 10} more")
        
        comp_choice = input("\nComparison strategy:\n   1) All pairwise\n   2) vs Control combination\nChoice (1/2): ")
        
        if comp_choice == '1':
            # All pairwise combinations
            for i in range(len(combo_options)):
                for j in range(i+1, len(combo_options)):
                    comparisons.append(f"{combo_options[i]['string']} vs {combo_options[j]['string']}")
        
        elif comp_choice == '2':
            # vs Control combination
            print("\nAvailable combinations:")
            for i, option in enumerate(combo_options):
                print(f"   {i+1:2d}) {option['string']}")
            
            control_idx = int(input("Select control combination (number): ")) - 1
            if 0 <= control_idx < len(combo_options):
                control_combo = combo_options[control_idx]
                
                for option in combo_options:
                    if option != control_combo:
                        comparisons.append(f"{option['string']} vs {control_combo['string']}")
    
    elif choice == '3':
        # Custom factor selection
        print("Select which factors to include in comparisons:")
        selected_factors = []
        for factor in all_factors:
            include = input(f"Include {factor}? (y/n): ")
            if include.lower() == 'y':
                selected_factors.append(factor)
        
        if selected_factors:
            print(f"\nUsing factors: {selected_factors}")
            combo_options = generate_combination_options(metadata_info, selected_factors)
            
            print(f"\nFound {len(combo_options)} possible combinations:")
            for i, option in enumerate(combo_options):
                print(f"   {i+1:2d}) {option['string']}")
            
            comp_choice = input("\nComparison strategy:\n   1) All pairwise\n   2) vs Control\nChoice (1/2): ")
            
            if comp_choice == '1':
                for i in range(len(combo_options)):
                    for j in range(i+1, len(combo_options)):
                        comparisons.append(f"{combo_options[i]['string']} vs {combo_options[j]['string']}")
            elif comp_choice == '2':
                control_idx = int(input("Select control combination (number): ")) - 1
                if 0 <= control_idx < len(combo_options):
                    control_combo = combo_options[control_idx]
                    for option in combo_options:
                        if option != control_combo:
                            comparisons.append(f"{option['string']} vs {control_combo['string']}")
    
    elif choice == '4':
        # All pairwise comparisons with all factors
        combo_options = generate_combination_options(metadata_info, all_factors)
        for i in range(len(combo_options)):
            for j in range(i+1, len(combo_options)):
                comparisons.append(f"{combo_options[i]['string']} vs {combo_options[j]['string']}")
        print(f"✅ Generated {len(comparisons)} pairwise comparisons")
    
    # Show final comparisons
    print(f"\n📋 Generated {len(comparisons)} comparisons:")
    for i, comp in enumerate(comparisons):
        print(f"   {i+1:2d}) {comp}")
    
    # Enhanced selection step - only for large numbers
    if len(comparisons) > 10:
        selected_comparisons = select_specific_comparisons(comparisons)
        return selected_comparisons
    else:
        return comparisons

def create_universal_config():
    """Main function to create universal configuration for selected dataset"""
    print("="*60)
    print(" MULTI-DATASET RNA-SEQ PIPELINE SETUP")
    print("="*60)
    
    # Step 1: Detect multiple datasets
    data_info = detect_multi_dataset_structure()
    if data_info['datasets_found'] == 0:
        print("❌ No valid datasets found")
        return None
    
    # Step 2: Select dataset
    selected_study, study_info = select_dataset(data_info)
    
    # Step 3: Select files within study
    selected_files = select_files_for_study(study_info)
    
    # Step 4: Analyze metadata structure
    metadata_info, metadata_df = analyze_metadata(selected_files['metadata'])
    if metadata_info is None:
        return None
    
    # Step 5: Column deletion (only if using original metadata)
    if selected_files['metadata_type'] == 'original':
        cleaned_metadata_info, cleaned_metadata_df, deleted_columns = interactive_column_deletion(metadata_df, metadata_info)
        
        # Save cleaned metadata if columns were deleted
        if deleted_columns:
            cleaned_metadata_file = save_cleaned_metadata(
                cleaned_metadata_df,
                selected_files['metadata'],
                deleted_columns,
                study_info['path']
            )
            selected_files['metadata'] = cleaned_metadata_file
    else:
        # Using pre-cleaned metadata
        cleaned_metadata_info = metadata_info
        cleaned_metadata_df = metadata_df
        deleted_columns = []
    
    # Step 6: Define experimental design
    exp_design = get_experimental_design(cleaned_metadata_info, cleaned_metadata_df)
    
    # Step 6.5: Validate sample matching
    sample_id_col = exp_design['sample_id_column']
    if not validate_sample_matching(selected_files['counts'], selected_files['metadata'], sample_id_col):
        print("❌ Sample mismatch detected. Please check your files.")
        return None
    
    # Step 7: Get analysis parameters
    analysis_params = get_analysis_parameters()
    
    # Step 8: Define comparisons (with enhanced selection)
    comparisons = define_comparisons(cleaned_metadata_info, exp_design)
    
    # Combine all configurations
    full_config = {
        'study_name': selected_study,
        'input_files': selected_files,
        'experimental_design': exp_design,
        'analysis_parameters': analysis_params,
        'comparisons': comparisons,
        'metadata_info': cleaned_metadata_info,
        'deleted_columns': deleted_columns if deleted_columns else [],
        'output_structure': 'by_contrast_csv_only'
    }
    
    # Step 9: Save configuration
    config_file = save_configuration(full_config, selected_study)
    
    # Step 10: Display summary
    display_configuration_summary(full_config, config_file)
    
    return full_config

def save_configuration(config, study_name):
    """Save the configuration to YAML file"""
    print("\n=== SAVING CONFIGURATION ===")
    
    # Create config directory if it doesn't exist
    config_dir = Path("conf")
    config_dir.mkdir(exist_ok=True)
    
    # Generate config filename with study name
    safe_study_name = "".join(c for c in study_name if c.isalnum() or c in ['_', '-'])
    config_filename = config_dir / f"analysis_config_{safe_study_name}.yaml"
    
    # Prepare config for saving (convert paths to strings, etc.)
    save_config = {}
    for key, value in config.items():
        if isinstance(value, dict):
            save_config[key] = {k: str(v) if isinstance(v, Path) else v for k, v in value.items()}
        elif isinstance(value, list):
            save_config[key] = [str(item) if isinstance(item, Path) else item for item in value]
        elif isinstance(value, Path):
            save_config[key] = str(value)
        else:
            save_config[key] = value
    
    # Save to YAML
    with open(config_filename, 'w') as f:
        yaml.safe_dump(save_config, f, default_flow_style=False, indent=2)
    
    print(f"✅ Configuration saved to: {config_filename}")
    
    # Also save a summary text file
    summary_file = config_dir / f"config_summary_{safe_study_name}.txt"
    
    with open(summary_file, 'w') as f:
        f.write(f"RNA-seq Analysis Configuration Summary\n")
        f.write(f"Generated on: {pd.Timestamp.now()}\n")
        f.write(f"Study: {study_name}\n")
        f.write(f"="*50 + "\n\n")
        
        f.write(f"Input Files:\n")
        f.write(f"  Count file: {config['input_files']['counts']}\n")
        f.write(f"  Metadata file: {config['input_files']['metadata']}\n")
        f.write(f"  Metadata type: {config['input_files']['metadata_type']}\n\n")
        
        f.write(f"Experimental Design:\n")
        f.write(f"  Sample ID column: {config['experimental_design']['sample_id_column']}\n")
        f.write(f"  Primary factor: {config['experimental_design']['primary_factor']}\n")
        f.write(f"  Additional factors: {config['experimental_design']['additional_factors']}\n")
        f.write(f"  Batch factor: {config['experimental_design']['batch_factor']}\n\n")
        
        f.write(f"Analysis Parameters:\n")
        f.write(f"  P-value threshold: {config['analysis_parameters']['pvalue_threshold']}\n")
        f.write(f"  Adjusted p-value threshold: {config['analysis_parameters']['padj_threshold']}\n")
        f.write(f"  Log fold-change threshold: {config['analysis_parameters']['logfc_threshold']}\n")
        f.write(f"  Minimum counts: {config['analysis_parameters']['min_counts']}\n\n")
        
        f.write(f"Comparisons ({len(config['comparisons'])}):\n")
        for i, comp in enumerate(config['comparisons'], 1):
            f.write(f"  {i:2d}) {comp}\n")
        
        if config['deleted_columns']:
            f.write(f"\nDeleted Columns ({len(config['deleted_columns'])}):\n")
            for col in config['deleted_columns']:
                f.write(f"  - {col}\n")
    
    print(f"✅ Configuration summary saved to: {summary_file}")
    return str(config_filename)

def display_configuration_summary(config, config_file):
    """Display final configuration summary"""
    print("\n" + "="*60)
    print(" CONFIGURATION COMPLETE")
    print("="*60)
    
    print(f"\n📊 Study: {config['study_name']}")
    print(f"📁 Configuration file: {config_file}")
    
    print(f"\n📂 Input Files:")
    print(f"   Count data: {Path(config['input_files']['counts']).name}")
    print(f"   Metadata: {Path(config['input_files']['metadata']).name}")
    print(f"   Metadata type: {config['input_files']['metadata_type']}")
    
    print(f"\n⚙️  Experimental Design:")
    print(f"   Sample ID: {config['experimental_design']['sample_id_column']}")
    print(f"   Primary factor: {config['experimental_design']['primary_factor']}")
    if config['experimental_design']['additional_factors']:
        print(f"   Additional factors: {', '.join(config['experimental_design']['additional_factors'])}")
    if config['experimental_design']['batch_factor']:
        print(f"   Batch factor: {config['experimental_design']['batch_factor']}")
    
    print(f"\n📊 Analysis Settings:")
    params = config['analysis_parameters']
    print(f"   P-value threshold: {params['pvalue_threshold']}")
    print(f"   Adjusted p-value threshold: {params['padj_threshold']}")
    print(f"   Log fold-change threshold: {params['logfc_threshold']}")
    print(f"   Minimum counts: {params['min_counts']}")
    
    print(f"\n🔬 Comparisons: {len(config['comparisons'])} total")
    if len(config['comparisons']) <= 10:
        for i, comp in enumerate(config['comparisons'], 1):
            print(f"   {i:2d}) {comp}")
    else:
        for i, comp in enumerate(config['comparisons'][:5], 1):
            print(f"   {i:2d}) {comp}")
        print(f"   ... and {len(config['comparisons']) - 5} more")
    
    if config['deleted_columns']:
        print(f"\n🗑️  Deleted columns: {len(config['deleted_columns'])}")
        if len(config['deleted_columns']) <= 5:
            for col in config['deleted_columns']:
                print(f"   - {col}")
        else:
            for col in config['deleted_columns'][:3]:
                print(f"   - {col}")
            print(f"   ... and {len(config['deleted_columns']) - 3} more")
    
    print(f"\n🎉 Setup complete! Ready to run the pipeline.")
    print(f"📋 Review configuration file: {config_file}")
    print(f"🚀 Run pipeline with: nextflow run main.nf --{config['study_name']}=true -profile standard")

def validate_configuration(config):
    """Validate the generated configuration"""
    errors = []
    warnings = []
    
    # Check required keys
    required_keys = ['study_name', 'input_files', 'experimental_design',
                     'analysis_parameters', 'comparisons']
    for key in required_keys:
        if key not in config:
            errors.append(f"Missing required configuration key: {key}")
    
    # Check input files exist
    if 'input_files' in config:
        for file_type, file_path in config['input_files'].items():
            if file_type in ['counts', 'metadata']:
                if not Path(file_path).exists():
                    errors.append(f"{file_type} file not found: {file_path}")
    
    # Check experimental design
    if 'experimental_design' in config:
        exp_design = config['experimental_design']
        if not exp_design.get('sample_id_column'):
            errors.append("Sample ID column not specified")
        if not exp_design.get('primary_factor'):
            errors.append("Primary factor not specified")
    
    # Check comparisons
    if 'comparisons' in config:
        if not config['comparisons']:
            warnings.append("No comparisons defined")
        elif len(config['comparisons']) > 50:
            warnings.append(f"Large number of comparisons ({len(config['comparisons'])})")
    
    return errors, warnings

def main():
    """Main entry point"""
    try:
        # Create configuration
        config = create_universal_config()
        if config is None:
            print("❌ Configuration setup failed")
            return 1
        
        # Validate configuration
        errors, warnings = validate_configuration(config)
        if errors:
            print(f"\n❌ Configuration errors found:")
            for error in errors:
                print(f"   - {error}")
            return 1
        
        if warnings:
            print(f"\n⚠️  Configuration warnings:")
            for warning in warnings:
                print(f"   - {warning}")
        
        print(f"\n🎉 Configuration setup completed successfully!")
        return 0
        
    except KeyboardInterrupt:
        print(f"\n\n⚠️  Setup interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    exit(main())