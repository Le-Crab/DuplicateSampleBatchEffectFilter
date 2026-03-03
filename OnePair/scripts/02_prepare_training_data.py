#!/usr/bin/env python3
"""
Script: 01_prepare_training_data.py
Purpose: Load and preprocess labeled training data for variant quality classification.
         Transforms annotated VCF data into ML-ready feature table.

Dependencies: pandas, numpy, scipy, pathlib
"""

import sys
from pathlib import Path
from typing import Tuple, Dict, List, Optional

import pandas as pd
from scipy.stats import zscore


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
    silent_read: bool = True
) -> Tuple[pd.DataFrame, Dict]:
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

    # Clean and z-score normalize depth/allele depth columns
    if 'AD_ALT' in df.columns and 'AD_REF' in df.columns and 'DP' in df.columns:
        print("Cleaning missing values in AD_REF and AD_ALT")
        omitted_af_merge = df['AD_REF'].isna().sum() + df['AD_ALT'].isna().sum()
        df = df.dropna(subset=['AD_REF', 'AD_ALT'])

        print("Applying z-score normalization to DP, AD_ALT, AD_REF")
        df['DP'] = zscore(df['DP'])
        df['AD_ALT'] = zscore(df['AD_ALT'])
        df['AD_REF'] = zscore(df['AD_REF'])

    dict_stats = {
        'processed_chunks': current_chunk,
        'omitted_while_reading': omitted_count,
        'omitted_after_merge': omitted_af_merge
    }

    return df, dict_stats


def prepare_training_data(
    fp_input: str,
    lst_selected_samples: List[str],
    fp_output: str
) -> None:
    """
    Load labeled data, filter to selected samples, and save preprocessed features.
    """
    
    print("\nSUBSTEP 1: LOADING LABELED DATA")
    

    # Define columns to keep
    lst_columns_to_keep = [
        'SNV_ID', 'GINeighbors', 'ReasontoDisc', 'RepeatMasker',
        'WindowMaskerSdust', 'SegDup', 'CpGisland', 'SimpleRepeat', 'Mappability',
        'GT', 'DP', 'AD_REF', 'AD_ALT', 'SB', 'GQ', 'VDiscordant'
    ]

    lst_numeric_cols = ['Mappability', 'DP', 'AD_REF', 'AD_ALT', 'SB', 'GQ', 'ReasontoDisc']
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

    # Remove low quality VConcordant (ReasontoDisc == 100)
    lenprior = len(df)
    df = df[df['ReasontoDisc'] != 100]
    print(f"Removed {lenprior - len(df):,} low quality VConcordant rows")
    print(f"Remaining rows: {len(df):,}")

    
    print("\nSUBSTEP2: FILTERING TO SELECTED SAMPLE PAIR")
    

    # Extract sample ID from SNV_ID (last part after underscore)
    df['SampleID'] = df['SNV_ID'].str.split('_').str[-1]
    print(f"Sample IDs to keep: {lst_selected_samples}")

    # Filter to selected samples
    df = df[df['SampleID'].isin(lst_selected_samples)]
    print(f"Rows after sample filtering: {len(df):,}")

    # Drop temporary SampleID column
    df = df.drop(['SampleID'], axis=1)

    
    print("\nSUBSTEP3: CREATING LABEL FROM VDISCORDANT")
    

    # Create label from VDiscordant
    df['label'] = df['VDiscordant']
    print(f"label distribution:\n{df['label'].value_counts()}")

    # Drop columns not needed for training
    df = df.drop(columns=['ReasontoDisc', 'SNV_ID', 'VDiscordant'])

    
    print("\nSUBSTEP4: SAVING PREPROCESSED DATA")
    

    # Create output directory if needed
    Path(fp_output).parent.mkdir(parents=True, exist_ok=True)

    # Save to CSV
    df.to_csv(fp_output, index=False)
    print(f"Saved preprocessed data to: {fp_output}")
    print(f"Final shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    print(f"Memory usage: {df.memory_usage(deep=True).sum() / 1024**2:.2f} MB")


def main() -> int:
    """Main execution function."""

    WD = Path(__file__).parent.parent

    # Input should come from script 00 output (labeled txt)
    INPUT_PATH = WD / 'data' / 'pg1kg30x6.bisnv.neighbor2all.annot2.samplewise.V.10M_labeled.txt'
	

    # Selected samples for one pair (same logic as script-style constants)
    SELECTED_SAMPLES = ['ERR194147', 'NA12878']

    # Output path
    OUTPUT_PATH = WD / 'data' / (INPUT_PATH.stem + '_mlready.csv')

    print(f"Input: {INPUT_PATH}")
    print(f"Output: {OUTPUT_PATH}")
    print(f"Selected samples: {SELECTED_SAMPLES}")

    # Validate input file exists
    if not INPUT_PATH.exists():
        print(f"[e] Input file not found: {INPUT_PATH}")
        return 1

    try:
        prepare_training_data(
            fp_input=str(INPUT_PATH),
            lst_selected_samples=SELECTED_SAMPLES,
            fp_output=str(OUTPUT_PATH)
        )
        print("Data preparation completed successfully")
        return 0

    except Exception as e:
        print(f"Error during data preparation: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
