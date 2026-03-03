#!/usr/bin/env python3
"""

Notes:
- This script was generated to handle large .txt inputs via chunked processing.
- Step 1 reading and preprocessing
- Step 2 Adding ncessary columns for labeling
- Step 3 Assignment of labels
- Step 4 Exporting the labeled dataset as txt file
"""
from pathlib import Path
import numpy as np
import pandas as pd

#Define the Advanced loop reading function 
def read_and_process_chunks(
    filepath,
    columns_to_keep,
    chunk_size=500000,
    category_cols=None,
    bool_cols=None,
    numeric_cols=None):
    """
    Reads a large CSV file in chunks, processes columns, and returns a cleaned DataFrame.
    """

    if category_cols is None:
        category_cols = []
    if bool_cols is None:
        bool_cols = []
    if numeric_cols is None:
        numeric_cols = []

    chunks = []
    current_chunk = 0

    for chunk in pd.read_csv(filepath, usecols=columns_to_keep, low_memory=False, sep='\t', chunksize=chunk_size):
        current_chunk += 1
        print('processed chunks: ', current_chunk)
        chunk = chunk[chunk['GQ']!= 0] # remove rows without GQ entry
        for col in numeric_cols:
            if col in chunk.columns:
                chunk[col] = pd.to_numeric(chunk[col], errors='coerce')
        for col in category_cols:
            if col in chunk.columns:
                chunk[col] = chunk[col].astype('category')
        for col in bool_cols:
            if col in chunk.columns:
                chunk[col] = chunk[col].astype('bool')   
                
        chunks.append(chunk)

    # If no chunks were read, return an empty DataFrame with the requested columns
    if not chunks:
        return pd.DataFrame(columns=columns_to_keep)

    df = pd.concat(chunks, ignore_index=True)

    # convert pair, sample to category
    if 'Pair' in df.columns:
        df['Pair'] = df['Pair'].astype('category')
    if 'Sample' in df.columns:
        df['Sample'] = df['Sample'].astype('category')

    return df



if __name__ == "__main__":
    WD = Path(__file__).parent.parent
    INPUT_PATH = WD / 'data' / 'pg1kg30x6.bisnv.neighbor2all.annot2.samplewise.V.10M.txt' # 10M subset
    # INPUT_PATH = WD / 'data' / 'pg1kg30x6.bisnv.neighbor2all.annot2.samplewise.sample.V.txt' # full size replace later
    SUFFIX = '_labeled'
    OUTPUT_PATH = WD / 'data' / (INPUT_PATH.stem + SUFFIX + '.txt')

    print(f"[i] Input: {INPUT_PATH}")
    print(f"[i] Output: {OUTPUT_PATH}")

    #================================== Step 1 ======================================+#
    # Reading the data with chunked reading loop
    columns_to_keep = [
        'CHROM', 'POS', 'REF', 'ALT', 'GINeighbors', 'GIDiscordant',  
        'RepeatMasker', 'WindowMaskerSdust', 'SegDup', 'SimpleRepeat', 'CpGisland',
        'Mappability', 'Sample', 'GT', 'DP', 'AD', 'SB', 'GQ'
    ]

    numeric_cols = ['POS', 'Mappability', 'DP', 'SB', 'GQ']
    category_cols = ['CHROM', 'REF', 'ALT', 'GINeighbors', 'GT','AD', 'Sample' ]  
    bool_cols = [ 'GIDiscordant', 'RepeatMasker', 'WindowMaskerSdust', 'SegDup', 'SimpleRepeat', 'CpGisland']

    df = read_and_process_chunks(
        filepath= INPUT_PATH, # train&Test dataset  extract
        columns_to_keep=columns_to_keep,
        chunk_size=2_000_000,
        category_cols=category_cols,
        bool_cols=bool_cols,
        numeric_cols=numeric_cols
    )

    # Cleaning the data
    # Removing anything but Complete Biallelic SNVs in GT
    allowed_gts = {'0/0', '0/1', '1/1'}  # only complete biallelic SNVs are allowed
    df = df[df['GT'].isin(allowed_gts)] # remove all but allowed GTs:

    df['GINeighbors'] = df['GINeighbors'].astype('category') # convert it after all combinations are present
    # Split AD column into AD_REF and AD_ALT by comma
    if 'AD' in df.columns:
        ad_split = df['AD'].str.split(',', expand=True)
        df['AD_REF'] = pd.to_numeric(ad_split[0], errors='coerce')
        df['AD_ALT'] = pd.to_numeric(ad_split[1], errors='coerce')
        # Remove rows where conversion failed (NA in AD_REF or AD_ALT)
        df = df.dropna(subset=['AD_REF', 'AD_ALT'])
        df = df.drop(columns=['AD'])

    df = df.dropna()
    print(f"Loaded {len(df):,} rows using chunked reading.")


    #========================================= Step 2 ======================================== #
    # Adding additional columns for assignment and preaparing to go into labeling step
    #  1. Create additional columns
    #  2. load sample allocation file and merge it with the main dataframe -> df_paired
    #  3. Validation of Discordant SNVs in paired samples by pivoting GTs in each pair and comparing them
    #  4. Transfer VDiscordant info back to df_paired  
    #  5. Remove unpaired SNVs from df_paired using a list of all_na_IDs

    # 2.1.  Create additional columns
    df['ID'] = df['CHROM'].astype(str) + '_' + df['POS'].astype(str) + '_' + df['REF'].astype(str)  + '_' + df['ALT'].astype(str) 
    df['SNV_ID']  = df['ID'] + '_' + df['Sample'].astype(str)  # Unique ID for each sample and variant
    df['ADDR'] = round( df['AD_ALT'] / (df['AD_REF'] + df['AD_ALT']),4)  # Allele Depth Ratio
    df = df.drop(columns= ['CHROM', 'POS', 'REF', 'ALT'  ])
    # Sus metrics
    df['susDP'] = np.where(df['DP'] <= 8, 1, 0)  # Suspect if DP < 8
    df['susGQ'] = np.where(df['GQ'] <= 20, 1, 0)  # Suspect if GQ < 20
    df['susSB'] = np.where(df['SB'] >= 10, 1, 0)  # Suspect if SB >= 10 
    df['susADDR'] = np.where((df['ADDR'] > 0) & (df['ADDR'] <= 0.2), 1, 0)  # Suspect if 0 < ADDR <= 0.2
    
    # 2.2. Paired configuration
    # Adding pair information to the dataset Samples by using external csv file.
    SampleAllocation = pd.read_csv(WD / 'data' / '1kgppg_paired_sample_allocation.csv', sep =',')
    SampleAllocation['Pair'] = SampleAllocation['Pair'].astype('category')
    #SampleAllocation = SampleAllocation[SampleAllocation['Pair'] != 'unpaired']  # Remove unpaired samples
    df_paired = df.merge(SampleAllocation, on='Sample', how='left') # join by Sample
    df_paired['Pair'] = df_paired['Pair'].astype('category')
    print(f"\n{len(df_paired):,} rows before removing singleton samples from dataset")
    #remove unpaired SNVs from df
    df_paired = df_paired[df_paired['Pair'] != 'unpaired'].copy()  # keep only paired samples
    print(f"{len(df_paired):,} rows after leaving only paired samples")


    # 2.3 # pivot comparison of GTs for all pairs
    # pivot comparison of GTs for all pairs
    all_mismatched_IDs_paired = [] # here VDiscordant IDs will be stored
    all_na_IDs= [] # here NA(unpaired IIDs with only one GY per pair) for further removal will be stored

    for pair in df_paired['Pair'].unique():
        if pair == 'unpaired':
            continue
        df_topivot = df_paired[df_paired['Pair'] == pair]
        
        # pivot_table with ID, Sample, GT
        table = df_topivot.pivot(index='ID', columns=['Sample'], values='GT')
        # get column names of resulting table
        compared_samples = table.columns
            
        # transfer only rows with NA in any of samples into a separate df
        table_na = table[table.isna().any(axis=1)]
        # extract unpaired (NA) from a pair
        na_IDs = table_na.index.tolist()

        # attach sample name to NA IDs
        na_IDs = [f"{id}_{compared_samples[0]}" for id in na_IDs]
        # append to the overall NA list
        all_na_IDs.extend(na_IDs)
            
        # 'GTmismatch' column to compare GT in pair, consider only allowedGTs when comparing
        allowedGTs = ['0/0', '0/1', '1/1']
        table['GTmismatch'] = (table[compared_samples[0]] != table[compared_samples[1]]) & (table[compared_samples[0]].isin(allowedGTs) & table[compared_samples[1]].isin(allowedGTs))
        
        # save IDs where GTmismatch == True
        mismatched_IDs = table[table['GTmismatch']].index.tolist()
        
        
        # append Sample name to mismatched IDs
        mismatched_IDs_paired = [f"{id}_{compared_samples[0]}" for id in mismatched_IDs] 
        
        # append to the overall list
        all_mismatched_IDs_paired.extend(mismatched_IDs_paired)
        
        # print(f"Mismatched IDs for {pair}: {len(mismatched_IDs)}")
        print(f"Mismatched GT pairs in {pair}: {len(mismatched_IDs_paired)}")
        # print(f"Example of mismatched ID paired: {mismatched_IDs_paired[:5]}")

    print(f"Total IDs to be VDiscordant: {len(all_mismatched_IDs_paired)}")
    print(f"Total number of unpaired IDs in 6 paired samples: {len(all_na_IDs)}")

    # 2.4. Transfer VDiscordant from all_mismatched_IDs_paired list to df_paired
    df_paired['VDiscordant'] = df_paired['SNV_ID'].isin(all_mismatched_IDs_paired)
    print(f"{df_paired['VDiscordant'].sum():,} rows with VDiscordant == True")

    # 2.5 Remove unpaired SNVs using a list of all_na_IDs
    print(f"{len(df_paired):,} rows abefore removing unpaired SNVs")
    df_paired = df_paired[~df_paired['SNV_ID'].isin(all_na_IDs)].copy()
    print(f"{len(df_paired):,} rows after removing unpaired SNVs")

    # --------------------- Concordant Correction #1 -----------------------#
    # experimental corrections that were devised during the analysis
    # Corrections reduce the false positive signal generated by Concordant SNVS with high susSum.
    # This particular correction reduces the impact of high SB in ALT homozygous Concordant SNVs
    # for rows with VDiscordant == False and GT == 1/1, set susSum to 0
    df_paired.loc[(df_paired['VDiscordant'] == False) & (df_paired['GT'] == '1/1'), 'susSB' ] = 0
    # --------------------- ------------------------------------------------#

    df_paired['susSum'] = df_paired['susDP'] + df_paired['susGQ'] + df_paired['susSB'] + df_paired['susADDR']

    #=================================== Step 3 =================================== #
    # Assignment of label
    # 1. Compute group max and assign labels
    df_paired['max_sus'] = df_paired.groupby(['ID','Pair'], observed=False)['susSum'].transform('max')
    df_paired['ReasontoDisc'] = 0
    # 2. Keep only rows with VDiscordant == True for labeling If susSum equals the group max, set ReasentoDisc = susSum; else leave 0
    df_paired.loc[
        (df_paired['VDiscordant'] == True) &
        (df_paired['susSum'] == df_paired['max_sus']),
        'ReasontoDisc'
    ] = df_paired['susSum']  


    # ------------------  Concordant Correction #2 ----------------------------------------#
    # Concordant Paired SNVs with susSum >0 are to be assigned ReasontoDisc =100
    # Concordant SNVs with susSum >0 are to be removed from training data to reduce confusion for the model 
    df_paired.loc[
        (df_paired['VDiscordant'] == False) &
        (df_paired['susSum'] > 0),
        'ReasontoDisc'
    ] = 100
    #--------------------------------------------------------------------------------------#

    print('\nlabels assigned, continuing to exporting')

    #====================================== Step 4 ======================================# 
    # 1. create df_assigned
    
    # 4.1 Creating df_paired_assigned and df_paired_assigned_humanread
    df_paired_assigned_humanread = df_paired.sort_values(by=['ID', 'Pair']).reset_index(drop=True) # create humanread
    df_paired_assigned = df_paired.drop([ 'ID', 'Sample', 'GIDiscordant', 'ADDR', 'susDP', 'susADDR', 'susSB', 'susGQ', 'max_sus', 'susSum', 'Pair'], axis =1)
    print('Columns in MLready dataset:')
    print(df_paired_assigned.columns.tolist())

 
    #Arrange by Sample and IDT
    # save to txt file
    df_paired_assigned.to_csv(OUTPUT_PATH, index=False, sep='\t')
    print('Labeling is done')
    print(f"Saved labeled data to: {OUTPUT_PATH}")
