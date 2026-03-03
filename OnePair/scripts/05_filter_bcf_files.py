#!/usr/bin/env python3
"""
Script: 05_filter_bcf_files.py
Purpose: Filter single-sample BCF/VCF.GZ files using predicted erroneous variant lists.
         Uses bcftools to exclude predicted low-quality variants from VCF files.

Dependencies: subprocess, pathlib
"""

import subprocess
import sys
from pathlib import Path
from typing import List, Tuple
# =============================================================================
# BCF FILTERING FUNCTIONS
# =============================================================================
def filter_single_bcf(
    fp_vcf: Path,
    fp_exclude_list: Path,
    dir_output: Path,
    suffix: str,
    bcftools: str = "bcftools") -> Tuple[bool, str]:
    """
    Filter a single VCF/BCF file by excluding variants in the provided list.
    
    Parameters
    ----------
    fp_vcf : Path
        Path to input VCF.GZ or BCF file
    fp_exclude_list : Path
        Path to file containing variant IDs to exclude (one per line)
    dir_output : Path
        Output directory for filtered BCF
    suffix : str
        Suffix to add to output filename
    bcftools : str
        Path to bcftools executable
    
    Returns
    -------
    Tuple[bool, str]
        Success status and output file path
    """
    sample_name = fp_vcf.stem
    if sample_name.endswith('.vcf'):
        sample_name = sample_name[:-4]  # Remove .vcf from .vcf.gz
    
    fp_output = dir_output / f"{sample_name}.{suffix}_Filtered.bcf"
    
    # Check if exclude list is empty
    if fp_exclude_list.stat().st_size == 0:
        print("  [w] Exclude list is empty. Copying original file.")
        # Copy original without filtering
        cmd = [
            bcftools, 'view',
            '-O', 'b',  # Output BCF format
            '-o', str(fp_output),
            str(fp_vcf)
        ]
    else:
        # Filter using bcftools
        cmd = [
            bcftools, 'view',
            '-e', f'ID=@{fp_exclude_list}',  # Exclude variants in file
            '-O', 'b',  # Output BCF format
            '-o', str(fp_output),
            str(fp_vcf)
        ]
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        return True, str(fp_output)
    except subprocess.CalledProcessError as e:
        msg = e.stderr.strip() if isinstance(e.stderr, str) else str(e)
        print(f"  [e] bcftools failed: {msg}")
        return False, ""


def index_bcf(fp_bcf: str, bcftools: str = "bcftools") -> bool:
    """
    Create index for BCF file using bcftools.
    
    Parameters
    ----------
    fp_bcf : str
        Path to BCF file
    bcftools : str
        Path to bcftools executable
    
    Returns
    -------
    bool
        Success status
    """
    cmd = [bcftools, 'index', '-f', fp_bcf]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError as e:
        msg = e.stderr.strip() if isinstance(e.stderr, str) else str(e)
        print(f"[e] Failed to index {fp_bcf}: {msg}")
        return False


def batch_filter_bcf_files(
    dir_vcf: str,
    dir_predictions: str,
    dir_output: str,
    suffix: str,
    vcf_pattern: str = "*.vcf.gz",
    bcftools: str = "bcftools") -> None:
    """
    Filter all VCF/BCF files in directory using corresponding prediction files.
    
    Parameters
    ----------
    dir_vcf : str
        Directory containing VCF.GZ or BCF files
    dir_predictions : str
        Directory containing prediction files (*_erroneous.txt)
    dir_output : str
        Output directory for filtered BCF files
    suffix : str
        Suffix for output files
    vcf_pattern : str
        File pattern to match (default: *.vcf.gz)
    bcftools : str
        Path to bcftools executable
    """
    print("\nSUBSTEP 1: DISCOVERING VCF FILES")
    
    dir_vcf = Path(dir_vcf)
    lst_vcf_files = sorted(dir_vcf.glob(vcf_pattern))
    print(f"Found {len(lst_vcf_files)} VCF files in {dir_vcf}")
    
    if not lst_vcf_files:
        print(f"[e] No files matching pattern '{vcf_pattern}' found")
        return
    
    dir_predictions = Path(dir_predictions)
    dir_output = Path(dir_output)
    dir_output.mkdir(parents=True, exist_ok=True)

    print("\nSUBSTEP 2: FILTERING VCF FILES")
    print(f"Output directory: {dir_output}")
    
    processed = 0
    skipped = 0
    
    for idx, fp_vcf in enumerate(lst_vcf_files, 1):
        sample_name = fp_vcf.stem
        if sample_name.endswith('.vcf'):
            sample_name = sample_name[:-4]
        
        if sample_name.startswith("sample_"):
            sample_name = sample_name.removeprefix("sample_")

        
        print(f"\n[{idx}/{len(lst_vcf_files)}] Processing: {fp_vcf.name}")
        
        # Find corresponding prediction file
        fp_pred = dir_predictions / f"{sample_name}_erroneous.txt"
        
        if not fp_pred.exists():
            print(f"  [w] Prediction file not found: {fp_pred.name}")
            skipped += 1
            continue
        
        # Count variants to exclude
        with open(fp_pred, 'r') as f:
            num_exclude = sum(1 for line in f if line.strip())
        
        print(f"  Variants to exclude: {num_exclude:,}")
        
        # Filter VCF
        success, fp_output = filter_single_bcf(
            fp_vcf=fp_vcf,
            fp_exclude_list=fp_pred,
            dir_output=dir_output,
            suffix=suffix,
            bcftools=bcftools
        )
        
        if success:
            print(f"  ✓ Created: {Path(fp_output).name}")
            
            # Index the filtered BCF
            print("  Indexing BCF...")
            if index_bcf(fp_output, bcftools=bcftools):
                print(f"  ✓ Indexed: {Path(fp_output).name}")
            else:
                print(f"  [w] Failed to index: {Path(fp_output).name}")
            
            processed += 1
        else:
            print(f"  [e] Failed to filter {fp_vcf.name}")
            skipped += 1
    
    print("\nSUBSTEP 3: FILTERING SUMMARY")
    print(f"Files processed: {processed}")
    print(f"Files skipped:   {skipped}")
    print(f"Total files:     {len(lst_vcf_files)}")


# =============================================================================
# MAIN
# =============================================================================
def main() -> int:
    """Main execution function."""
    WD = Path(__file__).parent.parent

    DIR_VCF = WD / 'data' / 'single_sample_vcfs'
    DIR_PREDICTIONS = WD / 'data' / 'erroneous_ids_single_sample_txts'
    DIR_OUTPUT = WD / 'data' / 'filtered_vcf' / 'filtered_single_sample_vcfs'

    SUFFIX = 'SB-OnePair-VDVC'
    VCF_PATTERN = '*.vcf.gz'
    BCFTOOLS = '/xcatopt/bcftools-1.18/bcftools'  # Adjust if bcftools is in PATH or located elsewhere

    print(f"VCF directory: {DIR_VCF}")
    print(f"Predictions directory: {DIR_PREDICTIONS}")
    print(f"Output directory: {DIR_OUTPUT}")
    print(f"Suffix: {SUFFIX}")
    print(f"VCF pattern: {VCF_PATTERN}")
    print(f"bcftools: {BCFTOOLS}")

    # Validate inputs
    if not DIR_VCF.exists():
        print(f"[e] VCF directory not found: {DIR_VCF}")
        return 1

    if not DIR_PREDICTIONS.exists():
        print(f"[e] Predictions directory not found: {DIR_PREDICTIONS}")
        return 1

    # Check bcftools availability
    try:
        subprocess.run([BCFTOOLS, '--version'], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"[e] bcftools not found or not executable: {BCFTOOLS}")
        return 1
    
    try:
        batch_filter_bcf_files(
            dir_vcf=str(DIR_VCF),
            dir_predictions=str(DIR_PREDICTIONS),
            dir_output=str(DIR_OUTPUT),
            suffix=SUFFIX,
            vcf_pattern=VCF_PATTERN,
            bcftools=BCFTOOLS
        )
        
        print("BCF filtering completed successfully")
        return 0
        
    except Exception as e:
        print(f"Error during BCF filtering: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
