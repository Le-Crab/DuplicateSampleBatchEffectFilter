#!/usr/bin/env python3
"""
Script: 02rel_prepare_training_data.py
Purpose: Load and preprocess labeled training data for variant quality classification.
         Transforms labeled txt into mlready txt


Dependencies: pandas, numpy, scipy, pathlib
"""

import sys
from pathlib import Path
from typing import Tuple, Dict, List, Optional

import pandas as pd
from scipy.stats import zscore

# I/O
WD = Path(__file__).parent.parent

# Input should come from script 00 output (labeled txt)
#INPUT_PATH = WD / 'data' / 'pgds1kgp.bisnv.neighbor2.annot.sw.labeled.subset1M.txt' 
INPUT_PATH = WD / 'data' / 'pgds1kgp.bisnv.neighbor2.annot.sw.labeled.txt'

# Output path
OUTPUT_PATH = WD / 'data' / (INPUT_PATH.stem + '.mlready.txt')
print(f"Input: {INPUT_PATH}")
print(f"Output: {OUTPUT_PATH}")


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================
def advanced_chunked_reading(
    fp_input: str,
    lst_usecols: List[str],
    lst_numeric_cols: Optional[List[str]] = None,
    lst_category_cols: Optional[List[str]] = None,
    lst_bool_cols: Optional[List[str]] = None,
    filter_notna_col: Optional[str] = None,
    chunk_size: int = 2_000_000,
    sep: str = '\t',
    silent_read: bool = True) -> Tuple[pd.DataFrame, Dict]:
    """
    Reads a large delimited file in chunks with dtype conversions and optional filtering.
    """
    lst_chunks = []
    omitted_count = 0
    current_chunk = 0
    omitted_af_merge = 0

    print(f"Reading file in chunks: {fp_input}")

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

        # Optional filter
        if filter_notna_col and filter_notna_col in chunk.columns:
            before = len(chunk)
            chunk = chunk[chunk[filter_notna_col].notna()]
            omitted_count += (before - len(chunk))

        if not silent_read:
            print(f"Processed chunk: {current_chunk}")

        lst_chunks.append(chunk)

    if not lst_chunks:
        print("[w] No data chunks read!")
        return pd.DataFrame(), {
            'processed_chunks': 0,
            'omitted_while_reading': 0,
            'omitted_after_merge': 0
        }

    df = pd.concat(lst_chunks, ignore_index=True)

    # Convert GINeighbors to category
    if 'GINeighbors' in df.columns:
        df['GINeighbors'] = df['GINeighbors'].astype('category')

    dict_stats = {
        'processed_chunks': current_chunk,
        'omitted_while_reading': omitted_count,
        'omitted_after_merge': omitted_af_merge
    }

    return df, dict_stats


def prepare_training_data(
    fp_input: str,
    fp_output: str) -> None:
    """
    Load labeled data, filter to selected samples, and save preprocessed features.
    """
    
    print("\nSUBSTEP 1: LOADING LABELED DATA")
    

    # Define columns to keep
    lst_columns_to_keep = [
        'SNV_ID', 'GINeighbors', 'ReasontoDisc', 'RepeatMasker',
        'WindowMaskerSdust', 'SegDup', 'CpGisland', 'SimpleRepeat', 'Mappability',
        'DP', 'GT', 'ADDR', 'SB', 'GQ', 'VDiscordant'
    ]

    lst_numeric_cols = ['Mappability', 'DP', 'ADDR', 'SB', 'GQ', 'ReasontoDisc']
    lst_category_cols = ['GT']
    lst_bool_cols = ['RepeatMasker', 'WindowMaskerSdust', 'SegDup', 'SimpleRepeat',
                     'CpGisland', 'VDiscordant']

    # Read data (silent_read defaults to True)
    df, dict_stats = advanced_chunked_reading(
        fp_input=fp_input,
        lst_usecols=lst_columns_to_keep,
        chunk_size=2_000_000,
        sep='\t',
        lst_numeric_cols=lst_numeric_cols,
        lst_category_cols=lst_category_cols,
        lst_bool_cols=lst_bool_cols
    )

    print(f"Total rows loaded: {len(df):,}")
    print(f"Reading statistics: {dict_stats}")

    # Report ReasontoDisc distribution
    print(f"ReasontoDisc distribution:\n{df['ReasontoDisc'].value_counts()}")
    
    #======================= CORRECTIONS AND FEATURE ENGINEERING ==============================#
    # Remove low quality VConcordant (ReasontoDisc == 100)
    lenprior = len(df)
    df = df[df['ReasontoDisc'] != 100]
    print(f"Removed {lenprior - len(df):,} low quality VConcordant rows")
    print(f"Remaining rows: {len(df):,}")
    
    print("\nSUBSTEP2: CREATING LABEL FROM VDISCORDANT")
    
    # Create label from VDiscordant
    df['Label'] = df['VDiscordant']
    print(f"label distribution:\n{df['Label'].value_counts()}")
    pos_label_ratio = df['Label'].value_counts().get(1, 0) / len(df)
    print(f"Positive label ratio: {pos_label_ratio:.4f}")

    # Z-score normalization for Depth
    if 'DP' in df.columns:
        df['DP'] = zscore(df['DP'])
    if 'AD_ALT' in df.columns:
        df['AD_ALT'] = zscore(df['AD_ALT'])
    

    # Drop columns not needed for training
    df = df.drop(columns=['ReasontoDisc', 'SNV_ID', 'VDiscordant'])

    # =========================================================================================#
    print("\nSUBSTEP3: SAVING PREPROCESSED DATA")
    

    # Create output directory if needed
    Path(fp_output).parent.mkdir(parents=True, exist_ok=True)

    # Save to CSV
    df.to_csv(fp_output, index=False, sep='\t')
    print(f"Saved preprocessed data to: {fp_output}")
    print(f"Final shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print(f"Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")


def main() -> int:
    """Main execution function."""

    # Validate input file exists
    if not INPUT_PATH.exists():
        print(f"[e] Input file not found: {INPUT_PATH}")
        return 1

    try:
        prepare_training_data(
            fp_input=str(INPUT_PATH),
            fp_output=str(OUTPUT_PATH)
        )
        print("Data preparation completed successfully")
        return 0

    except Exception as e:
        print(f"Error during data preparation: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
