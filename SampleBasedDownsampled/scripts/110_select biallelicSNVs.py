#!/usr/bin/env python3
"""
Select biallelic SNVs

Behavior:
- Finds files matching: Quartet_DNA_ILM_Nova_ARD_*.vcf.gz
- extracts bbiallelic SNVs with bcftools
-v snps - small variants only
-m2 - minimal number of alleles  = 2
-M2 - maximum number of alleles = 2
"""

from pathlib import Path
import glob
import subprocess
import re

bcftools = '/xcatopt/bcftools-1.18/bcftools'  # Adjust if bcftools is in PATH or located elsewhere
WD = Path(__file__).parent
DIR_VCF = WD /'data'  # input dir with VCF files
DIR_OUTPUT = DIR_VCF # output dir for biallelic SNVs


def main():
    # iterate over files matching the *.neighbor?all.annot.txt in WD
    
    pattern = "pgds1kgp.bcf" # pattern starts with 
    files = sorted(DIR_VCF.glob(pattern))
    print(f"Found files: {len(files)}")
    if not files:
        print(f"No files found matching pattern: {pattern}")
        return
    
    for file in files:
        print(f"Processing: {Path(file).name}")
        # Process each file
        cmd = [
            bcftools,
            "view",
            "-v",
            "snps",
            "-m2",
            "-M2",
            "-Oz",
            "-o",
            DIR_OUTPUT / str(file).replace('.vcf.gz', '.bial.vcf.gz'),
            file
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)



if __name__ == "__main__":
    main()
