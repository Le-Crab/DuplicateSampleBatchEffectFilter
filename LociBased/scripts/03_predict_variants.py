#!/usr/bin/env python3
"""
Script: 03_predict_variants.py
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
    silent_read: bool = True) -> Tuple[pd.DataFrame, Dict]:
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
        chunksize=chunk_size    ):
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
    df['GIMeanDepth'] = zscore(df['GIMeanDepth'])
    df['GIMeanDepthAlt'] = zscore(df['GIMeanDepthAlt'])
    

    dict_stats = {
        'processed_chunks': current_chunk,
        'omitted_after_merge': omitted_af_merge
    }
    
    return df, dict_stats


# =============================================================================
# PREDICTION FUNCTION
# =============================================================================
def predict_locibased(
    fp_input: Path,
    fp_model: Path,
    dir_output: Path) -> Tuple[int, int]:
    """
    Process multisample exported file and save predicted erroneous variant IDs.
    
    Parameters
    ----------
    fp_input : Path
        Path to input TXT file
    model_xgb : xgb.XGBClassifier
        Trained XGBoost model
    dir_output : Path
        Output directory for prediction files
    
    Returns
    -------
    Tuple[int, int]
        Total variants processed, number predicted as erroneous
    """
    
    
    # Define columns to load (ID column + features)
    lst_columns_to_keep =  [
    'POS', 'CHROM', 'REF', 'ALT', 
    'GIMeanDepth', 'GIMeanDepthAlt', 'GINeighbors', 'GIStrandBias', 'Mappability', 
    'RepeatMasker', 'WindowMaskerSdust', 'SegDup', 'SimpleRepeat', 'CpGisland']
    
    lst_numeric_cols = ['POS', 'GIMeanDepth', 'GIMeanDepthAlt', 'Mappability', 'GIStrandBias']
    lst_category_cols = ['CHROM', 'REF', 'ALT', 'GINeighbors']
    lst_bool_cols = [ 'RepeatMasker', 'WindowMaskerSdust', 'SegDup', 'SimpleRepeat', 'CpGisland']
    
    # Load sample data
    df, dict_stats = advanced_chunked_reading(
        fp_input=str(fp_input),
        lst_usecols=lst_columns_to_keep,
        chunk_size=2_000_000,
        sep='\t',
        lst_numeric_cols=lst_numeric_cols,
        lst_category_cols=lst_category_cols,
        lst_bool_cols=lst_bool_cols,
        silent_read=False
    )
    
    if df.empty:
        print(f"  [w] No data loaded for {fp_input}")
        return 0, 0
    
    total_variants = len(df)
    # Create ID column
    df['ID'] = df['CHROM'].astype(str) + '_' + df['POS'].astype(str) + "_" + df['REF'].astype(str) + '_' + df['ALT'].astype(str)
    # Separate IDs and features
    df_ids = df['ID']
    df_features = df.drop(['ID', 'CHROM', 'POS','REF', 'ALT'], axis=1)
    
    #load model
    model_xgb = xgb.XGBClassifier()
    model_xgb.load_model(fp_model)
    print(f"Loaded model from: {fp_model}")

    # Make predictions
    y_pred = model_xgb.predict(df_features)
    
    # Filter to erroneous predictions (positive class)
    # XGBoost predict() returns integer labels for binary classification
    mask_erroneous = (y_pred == 1)  
    lst_erroneous_ids = df_ids[mask_erroneous].tolist()
    num_erroneous = len(lst_erroneous_ids)
    
    # Save erroneous IDs to file with  "_erroneous_ids.txt" suffix
    fp_output = dir_output / f"{str(fp_input.stem)}_erroneous_ids.txt"
    Path(fp_output).parent.mkdir(parents=True, exist_ok=True)
    with open(fp_output, 'w') as f:
        for variant_id in lst_erroneous_ids:
            f.write(f"{variant_id}\n")
    
    return total_variants, num_erroneous


# =============================================================================
# MAIN
# =============================================================================
def main() -> int:
    """Main execution function."""
    WD = Path(__file__).parent.parent

    # Inputs/outputs
    # DIR_INPUT = WD / 'data' / 'single_sample_txts' / 'subset' # for dev; replace with 'single-sample-txts' for full run
    INPUT = WD / 'data' / 'yapobbj.bisnv.nocompneighbor2.annot2.simple.txt' # for full run
    MODEL_DIR = WD / 'models'
    DIR_OUTPUT = WD / 'data' / 'erroneous_ids'

    # Pick newest booster.json by default
    lst_models = sorted(MODEL_DIR.glob('*_booster.json'), key=lambda p: p.stat().st_mtime)
    FP_MODEL = lst_models[-1] if lst_models else None

    print(f"Input file: {INPUT}")
    print(f"Model directory: {MODEL_DIR}")
    print(f"Output directory: {DIR_OUTPUT}")
    # print(f"Will be saved to: {DIR_OUTPUT /  INPUT.stem / "_erroneous_ids.txt"}")
    # Validate inputs
    if not INPUT.exists():
        print(f"[e] Input file not found: {INPUT}")
        return 1

    if FP_MODEL is None or not FP_MODEL.exists():
        print(f"[e] No model found in: {MODEL_DIR} (expected '*_booster.json')")
        return 1
    
    try:
        predict_locibased(
            fp_input=INPUT,
            fp_model=str(FP_MODEL),
            dir_output=DIR_OUTPUT
        )
        
        print("Prediction completed successfully")
        return 0
        
    except Exception as e:
        print(f"Error during prediction: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
