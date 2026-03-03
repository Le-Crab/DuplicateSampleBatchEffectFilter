#!/usr/bin/env python3
"""
Script: 04_filter_bcf_file.py
Purpose: Filter ONE multi-sample VCF/BCF file using a list of predicted
         discordant variant IDs.

Dependencies: subprocess, pathlib
"""

import subprocess
import sys
from pathlib import Path
from typing import Tuple


def _check_bcftools(bcftools: str) -> bool:
    try:
        subprocess.run([bcftools, "--version"], check=True, capture_output=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


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


def filter_one_bcf(
    fp_in: Path,
    fp_discordant_ids: Path,
    fp_out: Path,
    bcftools: str,) -> Tuple[bool, str]:
    """Filter fp_in by excluding IDs in fp_discordant_ids."""
    fp_out.parent.mkdir(parents=True, exist_ok=True)

    # If exclude list is empty, still write an output for downstream steps.
    if fp_discordant_ids.stat().st_size == 0:
        print("[w] Discordant IDs list is empty; copying input via bcftools view")
        cmd = [bcftools, "view", "-Ob", "-o", str(fp_out), str(fp_in)]
    else:
        # Match the notebook expression style.
        cmd = [
            bcftools,
            "filter",
            "-e",
            f"ID =@{fp_discordant_ids}",
            str(fp_in),
            "-Ob",
            "-o",
            str(fp_out),
        ]

    try:
        subprocess.run(cmd, check=True)
        return True, str(fp_out)
    except subprocess.CalledProcessError as e:
        return False, f"bcftools failed: {e}"


# =============================================================================
# MAIN
# =============================================================================
def main() -> int:
    """Main execution function."""
    WD = Path(__file__).parent.parent

    fp_in = WD / "data" / "yapobbj.bisnv.nocompneighbor2.annot2.vcf.gz"
    fp_discordant_ids = WD / "data" / "erroneous_ids"/ "yapobbj.bisnv.nocompneighbor2.annot2.simple_erroneous_ids.txt"
    fp_out = WD / "data" / "filtered_vcf" / f"{fp_in.stem}_LociBased_filtered.bcf"
    
    #Toolpath
    bcftools = "/xcatopt/bcftools-1.18/bcftools" # 

    print(f"Input BCF: {fp_in}")
    print(f"Discordant IDs: {fp_discordant_ids}")
    print(f"Output BCF: {fp_out}")
    print(f"bcftools: {bcftools}")

    if not fp_in.exists():
        print(f"[e] Input BCF not found: {fp_in}")
        return 1

    if not fp_discordant_ids.exists():
        print(f"[e] Discordant IDs list not found: {fp_discordant_ids}")
        return 1

    if not _check_bcftools(bcftools):
        print(f"[e] bcftools not found or not executable: {bcftools}")
        return 1

    # Report size of exclude list
    with open(fp_discordant_ids) as f:
        n_ex = sum(1 for line in f if line.strip())
    print(f"IDs to exclude: {n_ex:,}")

    ok, msg = filter_one_bcf(fp_in, fp_discordant_ids, fp_out, bcftools=bcftools)
    if not ok:
        print(f"[e] {msg}")
        return 1

    print("Indexing output BCF...")
    if index_bcf(str(fp_out), bcftools=bcftools):
        print("Indexed")
    else:
        print("[w] Indexing failed")

    print("Filtering completed successfully")
    return 0


if __name__ == '__main__':
    sys.exit(main())
