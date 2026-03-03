#!/usr/bin/env bash

# =============================================================================
# DESCRIPTION
# =============================================================================
# Pipeline to transform BCF files into PCA, Manhattan, and QQ plots
# 
# Requirements:
#   1. BCF file (2 projects merged)
#   2. PHENO file (tab-delimited: FID IID PHENOTYPE)
#   3. PCA annotation files (one per sequencing center)
#
# Pipeline Steps:
#   0. Configuration: Set paths and parameters
#   1. Setup: Create directory structure
#   2. QC + PLINK: Convert BCF to PLINK format with QC filters
#   3. FlashPCA: Generate eigenvectors for PCA
#   4. GWAS: Association testing using PLINK 1.9
# =============================================================================

set -euo pipefail

# =============================================================================
# 0. CONFIGURATION - Paths and Parameters
# =============================================================================

# --- Base directories ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WD="$(cd "${SCRIPT_DIR}/.." && pwd)"   # repo root

# --- Input files ---
# By default, run the pipeline for BOTH a Control BCF and a Filtered BCF.
#


BCF_CONTROL="${WD}/data/BBJ1K+JAPON_122S_control.bcf" # Merged Control BCF
BCF_FILTERED="${WD}/data/filtered_vcf/SB-OnePair-VDVC_Merged.bcf" # Filtered merged BCF

# To process only ONE, comment out ONE of the variable definitions below.
# Default: both).
# BCF_INPUTS=()
# if [[ -n "${BCF_CONTROL-}" ]]; then
#   BCF_INPUTS+=("${BCF_CONTROL}")
# fi

if [[ -n "${BCF_FILTERED-}" ]]; then
  BCF_INPUTS+=("${BCF_FILTERED}")
fi

if [[ ${#BCF_INPUTS[@]} -eq 0 ]]; then
  echo "[e] No BCF inputs selected. Set BCF_CONTROL and/or BCF_FILTERED." >&2
  exit 1
fi

# Extra inputs from WD/data
PHENO_FILE="${WD}/data/phenotypes.txt"
REMOVER_FILE="${WD}/data/remover.txt"

# PCA annotation files
PCAANNOT_1KGPPG="${WD}/data/pcaannot_1kgppg.txt"
PCAANNOT_BBJ1KYAPP="${WD}/data/pcaanot_bbj1kyapp.txt"

# --- Output directories ---
OUT_DIR="${WD}/outputs"
SCRATCH_DIR="${OUT_DIR}/intermediate"
PLINK_SCRATCH="${SCRATCH_DIR}/plink"
PCA_SCRATCH="${SCRATCH_DIR}/pca"
PLOTS_DIR="${OUT_DIR}/plots"
GWAS_DIR="${OUT_DIR}/gwas"

# --- Tool paths ---
# Keep explicit paths for now (reproducibility), but allow override via env vars.
: "${PLINK2:=/xcatopt/plink2-20230621/plink2}"
: "${PLINK19:=/xcatopt/plink-1.90-beta6.21/plink}"
: "${FLASHPCA:=/xcatopt/flashpca/flashpca}"

# --- Verify critical inputs exist ---
echo "WD: ${WD}"
echo "BCF inputs:" \
  && for fp in "${BCF_INPUTS[@]}"; do echo "  - ${fp}"; done
echo "Pheno: ${PHENO_FILE}"
echo "Remover: ${REMOVER_FILE}"
echo "PCA annot 1KGPPG: ${PCAANNOT_1KGPPG}"
echo "PCA annot BBJ1KYAPP: ${PCAANNOT_BBJ1KYAPP}"

for fp in "${BCF_INPUTS[@]}" "${PHENO_FILE}" "${REMOVER_FILE}"; do
  if [[ ! -f "${fp}" ]]; then
    echo "[e] Required input not found: ${fp}" >&2
    exit 1
  fi
done

for tool in "${PLINK2}" "${PLINK19}" "${FLASHPCA}"; do
  if [[ ! -x "${tool}" ]]; then
    echo "[e] Required tool not executable: ${tool}" >&2
    exit 1
  fi
done

# =============================================================================
# 1. SETUP - Create Directory Structure
# =============================================================================

echo "SUBSTEP 1: SETUP DIRECTORIES"
mkdir -p "${PLINK_SCRATCH}" "${PCA_SCRATCH}" "${PLOTS_DIR}" "${GWAS_DIR}"
echo "Working directory: ${WD}"
echo "Directories prepared"

# =============================================================================
# 2-4. RUN PIPELINE FOR EACH INPUT BCF
# =============================================================================

for BCF_INPUT in "${BCF_INPUTS[@]}"; do
  OUTPREFIX="$(basename "${BCF_INPUT%.bcf}")"

  echo ""
  echo "SUBSTEP 2: PLINK QC + FORMAT CONVERSION (${OUTPREFIX})"
  "${PLINK2}" \
    --bcf "${BCF_INPUT}" \
    --make-bed \
    --pheno "${PHENO_FILE}" \
    --double-id \
    --allow-extra-chr \
    --not-chr 0,XY,X,Y,MT,decoyJRGv2 \
    --snps-only just-acgt \
    --vcf-half-call reference \
    --max-alleles 2 \
    --hwe 1e-6 \
    --geno-counts \
    --remove-fam "${REMOVER_FILE}" \
    --out "${PLINK_SCRATCH}/${OUTPREFIX}.raw" \
    2>&1 | tee "${PLINK_SCRATCH}/plink.${OUTPREFIX}.log"
  echo "PLINK processing finished (${OUTPREFIX})"

  echo "SUBSTEP 3: FLASH PCA (${OUTPREFIX})"
  "${FLASHPCA}" \
    --bfile "${PLINK_SCRATCH}/${OUTPREFIX}.raw" \
    --numthreads "$(nproc)" \
    --ndim 2 \
    --outpc "${PCA_SCRATCH}/${OUTPREFIX}.pc" \
    --outvec "${PCA_SCRATCH}/${OUTPREFIX}.vec" \
    --outval "${PCA_SCRATCH}/${OUTPREFIX}.val"
  echo "PCA eigenvectors generated (${OUTPREFIX})"

  echo "SUBSTEP 4: GWAS (PLINK 1.9) (${OUTPREFIX})"
  "${PLINK19}" \
    --bfile "${PLINK_SCRATCH}/${OUTPREFIX}.raw" \
    --assoc \
    --allow-no-sex \
    --pheno "${PHENO_FILE}" \
    --pheno-name PHENOTYPE \
    --out "${GWAS_DIR}/${OUTPREFIX}"
  echo "GWAS association test finished (${OUTPREFIX})"
done

echo ""
echo "SUBSTEP 5: DONE"
echo "Results summary:"
echo "  - PLINK files:  ${PLINK_SCRATCH}/"
echo "  - PCA outputs:  ${PCA_SCRATCH}/"
echo "  - GWAS results: ${GWAS_DIR}/"
echo "  - Plots:        ${PLOTS_DIR}/ (when generated)"
echo ""
