#!/usr/bin/env python3
"""
Script: 06_merge_bcf_files.py
Purpose: Merge filtered single-sample BCF files into one multi-sample BCF.
         Allows excluding corrupted/problematic samples from the merge.

Dependencies: subprocess, pathlib
"""

import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List
# =============================================================================
# BCF MERGING FUNCTIONS
# =============================================================================
def index_bcf(fp_bcf: Path, bcftools: str = "bcftools") -> bool:
    """
    Create index for BCF file using bcftools.
    
    Parameters
    ----------
    fp_bcf : Path
        Path to BCF file
    bcftools : str
        Path to bcftools executable
    
    Returns
    -------
    bool
        Success status
    """
    cmd = [bcftools, 'index', '-f', str(fp_bcf)]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        msg = e.stderr.strip() if isinstance(e.stderr, str) else str(e)
        print(f"[e] Failed to index {fp_bcf.name}: {msg}")
        return False


def merge_bcf_files(
    dir_bcf: str,
    fp_output: str,
    suffix: str,
    lst_exclude_samples: List[str] = None,
    bcftools: str = "bcftools",
    threads: int = 1
) -> None:
    """
    Merge filtered BCF files into a single multi-sample BCF.
    
    Parameters
    ----------
    dir_bcf : str
        Directory containing filtered BCF files
    fp_output : str
        Output path for merged BCF file
    suffix : str
        Suffix pattern to match BCF files (e.g., "SB-OnePair-VDVC")
    lst_exclude_samples : List[str], optional
        Sample IDs to exclude from merge
    bcftools : str
        Path to bcftools executable
    threads : int
        Number of threads for bcftools
    """
    print("\nSUBSTEP 1: DISCOVERING FILTERED BCF FILES")
    
    dir_bcf = Path(dir_bcf)
    lst_bcf_files = sorted(dir_bcf.glob("*.bcf"))
    print(f"Found {len(lst_bcf_files)} BCF files in {dir_bcf}")
    
    if not lst_bcf_files:
        print("[e] No BCF files found")
        return
    
    print("\nSUBSTEP 2: FILTERING BCF FILES BY SUFFIX AND EXCLUSIONS")
    
    # Filter BCFs by suffix
    lst_matching = []
    for fp in lst_bcf_files:
        if f".{suffix}_Filtered.bcf" in fp.name:
            lst_matching.append(fp)
    
    print(f"BCFs matching suffix '{suffix}': {len(lst_matching)}")
    
    # Apply sample exclusions
    if lst_exclude_samples is None:
        lst_exclude_samples = []
    
    lst_kept_files = []
    lst_omitted = []
    
    for fp in lst_matching:
        # Extract sample ID from filename
        # Format: {sample_id}.{suffix}.Filtered.bcf
        sample_id = fp.name.split(f".{suffix}.Filtered.bcf")[0]
        
        if sample_id in lst_exclude_samples:
            lst_omitted.append(sample_id)
            print(f"  Excluding: {sample_id}")
        else:
            lst_kept_files.append(fp)
    
    print(f"BCFs to keep:  {len(lst_kept_files)}")
    print(f"BCFs omitted:  {len(lst_omitted)}")
    if lst_omitted:
        print(f"Omitted samples: {', '.join(lst_omitted)}")
    
    if not lst_kept_files:
        print("[e] No BCF files remaining after filtering")
        return
    
    print("\nSUBSTEP 3: ENSURING ALL BCF FILES ARE INDEXED")
    
    for fp_bcf in lst_kept_files:
        fp_index = fp_bcf.with_suffix('.bcf.csi')
        if not fp_index.exists():
            print(f"Indexing: {fp_bcf.name}")
            if not index_bcf(fp_bcf, bcftools=bcftools):
                print(f"[e] Failed to index {fp_bcf.name}")
                return
        else:
            pass
    
    print("All BCF files are indexed")
    
    print("\nSUBSTEP 4: MERGING BCF FILES")
    
    # Create temporary file list for bcftools merge
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as fp_temp:
        for fp_bcf in lst_kept_files:
            fp_temp.write(f"{fp_bcf}\n")
        fp_filelist = fp_temp.name
    
    print(f"Created file list: {fp_filelist}")
    print(f"Merging {len(lst_kept_files)} BCF files...")
    print(f"Output: {fp_output}")
    
    # Merge using bcftools
    cmd = [
        bcftools, 'merge',
        '--file-list', fp_filelist,
        '--output', fp_output,
        '--output-type', 'b',  # BCF format
        '--threads', str(threads)
    ]
    
    try:
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        print("✓ Merge completed successfully")
        
        # Clean up temporary file list
        Path(fp_filelist).unlink()
        
        # Index the merged BCF
        print("\nSUBSTEP 5: INDEXING MERGED BCF")
        
        if index_bcf(Path(fp_output), bcftools=bcftools):
            print(f"✓ Indexed merged BCF: {fp_output}")
        else:
            print(f"[w] Failed to index merged BCF: {fp_output}")
        
        # Report statistics
        print("\nSUBSTEP 6: MERGE SUMMARY")
        print(f"Input BCFs:      {len(lst_kept_files)}")
        print(f"Excluded BCFs:   {len(lst_omitted)}")
        print(f"Output file:     {fp_output}")
        
        # Count samples in merged file
        cmd_count = [bcftools, 'query', '-l', fp_output]
        result = subprocess.run(cmd_count, capture_output=True, text=True, check=True)
        num_samples = len(result.stdout.strip().split('\n'))
        print(f"Samples in merged BCF: {num_samples}")
        
    except subprocess.CalledProcessError as e:
        msg = e.stderr.strip() if isinstance(e.stderr, str) else str(e)
        print(f"[e] bcftools merge failed: {msg}")
        Path(fp_filelist).unlink()
        raise


# =============================================================================
# MAIN
# =============================================================================
def main() -> int:
    """Main execution function."""
    WD = Path(__file__).parent.parent

    DIR_INPUT = WD / 'data' / 'filtered_vcf' / 'filtered_single_sample_vcfs'
    DIR_OUTPUT = WD / 'data' / 'filtered_vcf'

    SUFFIX = 'SB-OnePair-VDVC'
    FP_OUTPUT = DIR_OUTPUT / f"{SUFFIX}_Merged.bcf"

    # Optional config
    EXCLUDE_SAMPLES: List[str] = []
    BCFTOOLS = '/xcatopt/bcftools-1.18/bcftools'
    THREADS = 40

    print(f"Input directory: {DIR_INPUT}")
    print(f"Output file: {FP_OUTPUT}")
    print(f"Suffix: {SUFFIX}")
    if EXCLUDE_SAMPLES:
        print(f"Excluded samples: {', '.join(EXCLUDE_SAMPLES)}")
    print(f"bcftools: {BCFTOOLS}")
    print(f"Threads: {THREADS}")

    # Validate inputs
    if not DIR_INPUT.exists():
        print(f"[e] Input directory not found: {DIR_INPUT}")
        return 1

    # Check bcftools availability
    try:
        subprocess.run([BCFTOOLS, '--version'], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"[e] bcftools not found or not executable: {BCFTOOLS}")
        return 1

    # Ensure output directory exists
    FP_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        merge_bcf_files(
            dir_bcf=str(DIR_INPUT),
            fp_output=str(FP_OUTPUT),
            suffix=SUFFIX,
            lst_exclude_samples=EXCLUDE_SAMPLES,
            bcftools=BCFTOOLS,
            threads=THREADS
        )
        
        print("BCF merging completed successfully")
        return 0
        
    except Exception as e:
        print(f"Error during BCF merging: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
