"""
utils_pluto.py
==============
Handles:
- Pluto API downloads
- Local CSV loading
- Metadata factor detection
- Contrast generation
Follows Nextflow training best practices
"""
import os
import pandas as pd
import requests
from itertools import combinations
try:
    from plutobio import PlutoClient
    PLUTOBIO_AVAILABLE = True
except ImportError:
    PLUTOBIO_AVAILABLE = False
    print("Warning: plutobio not installed. Using direct API calls.")

# ------------------------------------------------------------
# Token Management
# ------------------------------------------------------------
def get_api_token():
    """Get API token from file or user input with validation"""
    token_file = "pluto_token.txt"
    
    # Check if token file exists
    if os.path.exists(token_file):
        with open(token_file, 'r') as f:
            stored_token = f.read().strip()
        if stored_token and stored_token != 'n':  # Skip invalid 'n' token
            print(f"✓ Found existing API token: {stored_token[:8]}...{stored_token[-8:]}")
            use_existing = input("Use existing token? (y/n): ").strip().lower()
            if use_existing in ('y', 'yes', ''):
                return stored_token
    
    # Get new token from user
    print("\nPlease enter a valid Pluto API token:")
    print("(You can get your token from https://pluto.bio/settings/api-tokens)")
    new_token = input("Token: ").strip()
    
    if new_token and new_token.lower() != 'n':
        # Save token for future use
        with open(token_file, 'w') as f:
            f.write(new_token)
        print("✓ Token saved for future use")
        return new_token
    else:
        raise RuntimeError("✗ No valid API token provided")

# ------------------------------------------------------------
# Load Data (Pluto or Local)
# ------------------------------------------------------------
def load_local_data(counts_path, metadata_path):
    counts = pd.read_csv(counts_path)
    metadata = pd.read_csv(metadata_path)
    return counts, metadata

def load_pluto_data(experiment_id, api_token=None):
    if not api_token:
        api_token = os.getenv("PLUTO_API_TOKEN")
    if not api_token:
        raise RuntimeError("✗ API token not provided")
    
    print(f"✓ Loading data for experiment: {experiment_id}")
    
    # First validate the token with a simple test call
    test_headers = {"Authorization": f"Token {api_token}"}
    test_url = "https://api.pluto.bio/lab/experiments/"
    
    try:
        test_response = requests.get(test_url, headers=test_headers)
        if test_response.status_code == 401:
            print("✗ API token is invalid or expired")
            print("Please get a new token from https://pluto.bio/settings/api-tokens")
            # Remove invalid token
            if os.path.exists("pluto_token.txt"):
                os.remove("pluto_token.txt")
            raise RuntimeError("Authentication failed - invalid API token")
    except requests.exceptions.RequestException as e:
        print(f"⚠ Could not validate token (network issue): {e}")
        print("Proceeding with data loading...")
    
    if PLUTOBIO_AVAILABLE:
        # Use PlutoClient SDK
        print("✓ Using PlutoClient SDK")
        client = PlutoClient(token=api_token)
        
        try:
            # This is a placeholder - actual PlutoClient methods may differ
            # You might need to adjust based on PlutoClient documentation
            sample_data = client.get_sample_data(experiment_id)
            assay_data = client.get_assay_data(experiment_id)
            
            metadata = pd.DataFrame(sample_data)
            counts = pd.DataFrame(assay_data)
            
        except Exception as e:
            print(f"⚠ PlutoClient failed, falling back to direct API: {e}")
            return _load_via_direct_api(experiment_id, api_token)
    else:
        # Use direct API calls
        return _load_via_direct_api(experiment_id, api_token)
    
    print(f"✓ Sample metadata shape: {metadata.shape}")
    print(f"✓ Assay counts shape: {counts.shape}")
    
    # Save to files in experiment ID folder
    data_folder = f"data/{experiment_id}"
    os.makedirs(data_folder, exist_ok=True)
    counts_path = f"{data_folder}/counts.csv"
    metadata_path = f"{data_folder}/metadata.csv"
    
    counts.to_csv(counts_path, index=False)
    metadata.to_csv(metadata_path, index=False)
    
    print(f"✓ Data saved:")
    print(f"  - Counts: {counts_path}")
    print(f"  - Metadata: {metadata_path}")
    
    return counts, metadata, counts_path, metadata_path

def _load_via_direct_api(experiment_id, api_token):
    """Load data using direct API calls based on Pluto documentation"""
    print("✓ Using direct API calls")
    
    # Correct headers format from Pluto documentation
    headers = {"Authorization": f"Token {api_token}"}
    limit = "10000"
    
    # Correct endpoint URLs from Pluto documentation
    assay_endpoint = f"https://api.pluto.bio/lab/experiments/{experiment_id}/assay-data/?limit={limit}"
    sample_endpoint = f"https://api.pluto.bio/lab/experiments/{experiment_id}/sample-data/?limit={limit}"
    
    try:
        print(f"✓ Fetching assay data...")
        assay_resp = requests.get(assay_endpoint, headers=headers)
        
        if assay_resp.status_code == 401:
            raise RuntimeError(f"✗ Authentication failed: Invalid API token for experiment {experiment_id}")
        elif assay_resp.status_code == 404:
            raise RuntimeError(f"✗ Experiment {experiment_id} not found or not accessible with this token")
        
        assay_resp.raise_for_status()
        assay_data = assay_resp.json()
        
        print(f"✓ Fetching sample data...")
        sample_resp = requests.get(sample_endpoint, headers=headers)
        
        if sample_resp.status_code == 401:
            raise RuntimeError(f"✗ Authentication failed: Invalid API token for experiment {experiment_id}")
        elif sample_resp.status_code == 404:
            raise RuntimeError(f"✗ Sample data for experiment {experiment_id} not found or not accessible")
        
        sample_resp.raise_for_status()
        sample_data = sample_resp.json()
        
    except requests.exceptions.HTTPError as e:
        if "401" in str(e):
            print("✗ Authentication failed - your API token is invalid or expired")
            print("Please get a new token from https://pluto.bio/settings/api-tokens")
            # Remove invalid token
            if os.path.exists("pluto_token.txt"):
                os.remove("pluto_token.txt")
        raise RuntimeError(f"API request failed: {e}")
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Network error: {e}")
    
    # Based on Pluto documentation, data is in "items" key
    if isinstance(assay_data, dict) and 'items' in assay_data:
        counts = pd.DataFrame(assay_data['items'])
    else:
        counts = pd.DataFrame(assay_data)
        
    if isinstance(sample_data, dict) and 'items' in sample_data:
        metadata = pd.DataFrame(sample_data['items'])
    else:
        metadata = pd.DataFrame(sample_data)
    
    print(f"✓ Sample metadata shape: {metadata.shape}")
    print(f"✓ Assay counts shape: {counts.shape}")
    
    # Save to files in experiment ID folder
    data_folder = f"data/{experiment_id}"
    os.makedirs(data_folder, exist_ok=True)
    counts_path = f"{data_folder}/counts.csv"
    metadata_path = f"{data_folder}/metadata.csv"
    
    counts.to_csv(counts_path, index=False)
    metadata.to_csv(metadata_path, index=False)
    
    print(f"✓ Data saved:")
    print(f"  - Counts: {counts_path}")
    print(f"  - Metadata: {metadata_path}")
    
    return counts, metadata, counts_path, metadata_path

# ------------------------------------------------------------
# Metadata Factor Detection
# ------------------------------------------------------------
def detect_factors(metadata_df):
    """Identify categorical fields (≤ 10 unique values)."""
    factors = []
    ignore_cols = {"sample", "sample_id", "id"} 
    
    for col in metadata_df.columns:
        unique_count = metadata_df[col].nunique()
        if col.lower() in ignore_cols:
            continue
        if unique_count <= 10:
            factors.append(col)
    
    if not factors:
        raise RuntimeError("✗ No categorical factors detected in metadata.")
    return factors

# ------------------------------------------------------------
# Contrast Generation
# ------------------------------------------------------------
def generate_comparisons(metadata_df, primary_factor):
    levels = metadata_df[primary_factor].unique()
    comps = []
    for a, b in combinations(levels, 2):
        comps.append(f"{primary_factor}:{a} vs {primary_factor}:{b}")
    return comps