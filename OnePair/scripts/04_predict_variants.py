#!/usr/bin/env python3
"""
Script: 04_predict_variants.py
Purpose: Apply trained XGBoost model to predict erroneous variants in single-sample files.
         Processes multiple samples in batch and outputs per-sample prediction lists.

Dependencies: pandas, numpy, xgboost, scipy
"""

import sys
from pathlib import Path
from typing import List, Tuple, Dict

import pandas as pd
import numpy as np
import xgboost as xgb
from scipy.stats import zscore
# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================
def advanced_chunked_reading(
    fp_input: str,
    lst_usecols: List[str] = None,
    chunk_size: int = 2_000_000,
    sep: str = '\t',
    lst_numeric_cols: List[str] = None,
    lst_category_cols: List[str] = None,
    lst_bool_cols: List[str] = None,
    silent_read: bool = True
) -> Tuple[pd.DataFrame, Dict]:
    """
    Read large delimited file in chunks with dtype conversions.
    
    Parameters
    ----------
    fp_input : str
        Path to input file
    lst_usecols : List[str], optional
        Columns to read
    chunk_size : int
        Number of rows per chunk
    sep : str
        Column separator
    lst_numeric_cols : List[str], optional
        Columns to convert to numeric
    lst_category_cols : List[str], optional
        Columns to convert to category
    lst_bool_cols : List[str], optional
        Columns to convert to boolean
    silent_read : bool
        Suppress chunk progress messages
    
    Returns
    -------
    Tuple[pd.DataFrame, Dict]
        DataFrame and statistics dictionary
    """
    lst_chunks = []
    current_chunk = 0
    omitted_af_merge = 0

    for chunk in pd.read_csv(
        fp_input, 
        usecols=lst_usecols, 
        low_memory=False, 
        sep=sep, 
        chunksize=chunk_size
    ):
        current_chunk += 1
        
        # Numeric conversions
        if lst_numeric_cols:
            for col in lst_numeric_cols:
                if col in chunk.columns:
                    chunk[col] = pd.to_numeric(chunk[col], errors='coerce')

        # Category conversions
        if lst_category_cols:
            for col in lst_category_cols:
                if col in chunk.columns:
                    chunk[col] = chunk[col].astype('category')

        # Boolean conversions
        if lst_bool_cols:
            for col in lst_bool_cols:
                if col in chunk.columns:
                    chunk[col] = chunk[col].astype('bool')
        
        if not silent_read:
            print(f"Processed chunk: {current_chunk}")

        lst_chunks.append(chunk)

    if not lst_chunks:
        return pd.DataFrame(), {'processed_chunks': 0, 'omitted_after_merge': 0}

    df = pd.concat(lst_chunks, ignore_index=True)
    
    # Convert GINeighbors to category
    if 'GINeighbors' in df.columns:
        df['GINeighbors'] = df['GINeighbors'].astype('category')

    # Clean and z-score normalize depth/allele depth columns
    if 'AD_ALT' in df.columns and 'AD_REF' in df.columns and 'DP' in df.columns:
        omitted_af_merge = df['AD_REF'].isna().sum() + df['AD_ALT'].isna().sum()
        df = df.dropna(subset=['AD_REF', 'AD_ALT'])
        
        df['DP'] = zscore(df['DP'])
        df['AD_ALT'] = zscore(df['AD_ALT'])
        df['AD_REF'] = zscore(df['AD_REF'])

    dict_stats = {
        'processed_chunks': current_chunk,
        'omitted_after_merge': omitted_af_merge
    }
    
    return df, dict_stats


# =============================================================================
# PREDICTION FUNCTION
# =============================================================================
def predict_single_sample(
    fp_sample: Path,
    model_xgb: xgb.XGBClassifier,
    lst_feature_names: List[str],
    dir_output: Path
) -> Tuple[int, int]:
    """
    Process single sample file and save predicted erroneous variant IDs.
    
    Parameters
    ----------
    fp_sample : Path
        Path to sample TXT file subjected to filtering
    model_xgb : xgb.XGBClassifier
        Trained XGBoost model
    lst_feature_names : List[str]
        Feature column names expected by model
    dir_output : Path
        Output directory for prediction files
    
    Returns
    -------
    Tuple[int, int]
        Total variants processed, number predicted as erroneous
    """
    sample_name = fp_sample.stem
    
    # Define columns to load (ID column + features)
    lst_columns_to_keep = ['ID'] + lst_feature_names
    
    lst_numeric_cols = ['Mappability', 'DP', 'AD_REF', 'AD_ALT', 'SB', 'GQ']
    lst_category_cols = ['GINeighbors', 'GT']
    lst_bool_cols = ['RepeatMasker', 'WindowMaskerSdust', 'SegDup', 'SimpleRepeat', 'CpGisland']
    
    # Load sample data
    df_sample, dict_stats = advanced_chunked_reading(
        fp_input=str(fp_sample),
        lst_usecols=lst_columns_to_keep,
        chunk_size=2_000_000,
        sep='\t',
        lst_numeric_cols=lst_numeric_cols,
        lst_category_cols=lst_category_cols,
        lst_bool_cols=lst_bool_cols,
        silent_read=True
    )
    
    if df_sample.empty:
        print(f"  [w] No data loaded for {sample_name}")
        return 0, 0
    
    total_variants = len(df_sample)
    
    # Separate IDs and features
    df_ids = df_sample['ID']
    df_features = df_sample[lst_feature_names]
    
    # Make predictions
    y_pred = model_xgb.predict(df_features)
    
    # Filter to erroneous predictions (positive class)
    # XGBoost predict() returns integer labels for binary classification
    mask_erroneous = (y_pred == 1)
    
    lst_erroneous_ids = df_ids[mask_erroneous].tolist()
    num_erroneous = len(lst_erroneous_ids)
    
    # Save erroneous IDs to file
    fp_output = dir_output / f"{sample_name}_erroneous.txt"
    with open(fp_output, 'w') as f:
        for variant_id in lst_erroneous_ids:
            f.write(f"{variant_id}\n")
    
    return total_variants, num_erroneous


def batch_predict_samples(
    dir_input: str,
    fp_model: str,
    dir_output: str,
    pattern: str = "*.txt"
) -> None:
    """
    Apply trained model to all sample files in input directory.
    
    Parameters
    ----------
    dir_input : str
        Directory containing sample TXT files
    fp_model : str
        Path to trained model JSON file
    dir_output : str
        Output directory for prediction files
    pattern : str
        File pattern to match (default: *.txt)
    """
    print("\nSUBSTEP 1: LOADING TRAINED MODEL")
    
    # Load model
    model_xgb = xgb.XGBClassifier()
    model_xgb.load_model(fp_model)
    print(f"Loaded model from: {fp_model}")
    
    # Get feature names from model
    lst_feature_names = model_xgb.get_booster().feature_names
    print(f"Model expects {len(lst_feature_names)} features")
    print(f"Features: {lst_feature_names}")
    
    print("\nSUBSTEP 2: DISCOVERING SAMPLE FILES")
    
    dir_input = Path(dir_input)
    lst_sample_files = sorted(dir_input.glob(pattern))
    print(f"Found {len(lst_sample_files)} sample files in {dir_input}")
    
    if not lst_sample_files:
        print(f"[e] No files matching pattern '{pattern}' found in {dir_input}")
        return
    
    print("\nSUBSTEP 3: PROCESSING SAMPLES")
    
    # Create output directory
    dir_output = Path(dir_output)
    dir_output.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {dir_output}")
    
    # Process each sample
    total_samples_processed = 0
    total_variants_processed = 0
    total_erroneous_predicted = 0
    
    for idx, fp_sample in enumerate(lst_sample_files, 1):
        print(f"\n[{idx}/{len(lst_sample_files)}] Processing: {fp_sample.name}")
        
        try:
            num_variants, num_erroneous = predict_single_sample(
                fp_sample=fp_sample,
                model_xgb=model_xgb,
                lst_feature_names=lst_feature_names,
                dir_output=dir_output
            )
            
            total_samples_processed += 1
            total_variants_processed += num_variants
            total_erroneous_predicted += num_erroneous
            
            if num_variants > 0:
                pct_erroneous = (num_erroneous / num_variants) * 100
                print(f"  Total variants: {num_variants:,}")
                print(f"  Erroneous:      {num_erroneous:,} ({pct_erroneous:.2f}%)")
            
        except Exception as e:
            print(f"  [e] Error processing {fp_sample.name}: {e}")
            continue
    
    print("\nSUBSTEP 4: PREDICTION SUMMARY")
    print(f"Samples processed:     {total_samples_processed}")
    print(f"Total variants:        {total_variants_processed:,}")
    print(f"Erroneous predicted:   {total_erroneous_predicted:,}")
    if total_variants_processed > 0:
        overall_pct = (total_erroneous_predicted / total_variants_processed) * 100
        print(f"Overall error rate:    {overall_pct:.2f}%")


# =============================================================================
# MAIN
# =============================================================================
def main() -> int:
    """Main execution function."""
    WD = Path(__file__).parent.parent

    # Inputs/outputs
    # DIR_INPUT = WD / 'data' / 'single_sample_txts' / 'subset' # for dev; replace with 'single-sample-txts' for full run
    DIR_INPUT = WD / 'data' / 'single_sample_txts' # for full run
    MODEL_DIR = WD / 'models'
    DIR_OUTPUT = WD / 'data' / 'erroneous_ids_single_sample_txts'
    PATTERN = '*.txt'

    # Pick newest booster.json by default
    lst_models = sorted(MODEL_DIR.glob('*_booster.json'), key=lambda p: p.stat().st_mtime)
    FP_MODEL = lst_models[-1] if lst_models else None

    print(f"Input directory: {DIR_INPUT}")
    print(f"Model directory: {MODEL_DIR}")
    print(f"Output directory: {DIR_OUTPUT}")
    print(f"Pattern: {PATTERN}")

    # Validate inputs
    if not DIR_INPUT.exists():
        print(f"[e] Input directory not found: {DIR_INPUT}")
        return 1

    if FP_MODEL is None or not FP_MODEL.exists():
        print(f"[e] No model found in: {MODEL_DIR} (expected '*_booster.json')")
        return 1
    
    try:
        batch_predict_samples(
            dir_input=str(DIR_INPUT),
            fp_model=str(FP_MODEL),
            dir_output=str(DIR_OUTPUT),
            pattern=PATTERN
        )
        
        print("Prediction completed successfully")
        return 0
        
    except Exception as e:
        print(f"Error during prediction: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
