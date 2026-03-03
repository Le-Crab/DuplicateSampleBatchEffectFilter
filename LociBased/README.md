# LociBased

This repository contains a simplified **loci-based** workflow featuring multissmple vcf file exports processing:

1. Preparing an ML-ready feature table from annotated vcf export
2. Training an XGBoost classifier
3. Predicting discordant/erroneous variants
4. Filtering a multi-sample VCF/BCF using the predicted IDs
5. Running PLINK + FlashPCA + GWAS, then making PCA / Manhattan / QQ plots

## Runnning the pipeline scripts

All steps are under `scripts/` and can be run sequentially via:
- `scripts/00_orchestration.py`

Or individually starting from the first:

- `scripts/01_prepare_training_data.py`
- `scripts/02_train_xgb_classifier.py`
- `scripts/03_predict_variants.py`
- `scripts/04_filter_bcf_file.py`
- `scripts/05_bcf2pca.sh`
- `scripts/06_vec2pca.R`
- `scripts/07_mhqq.R`

## Python requirements

The Python modules imported by the workflow scripts (`scripts/00`–`04`) are:

- `pandas`
- `numpy`
- `scipy`
- `xgboost`
- `scikit-learn`
- `matplotlib`

See also: `requirements.txt`.

## R requirements

The plotting scripts use:

- `ggplot2`
- `ggrepel`
- `dplyr`
- `CMplot`

## External CLI requirements

The downstream PCA/GWAS step requires these executables (paths configurable via env vars in `scripts/05_bcf2pca.sh`):

- `plink2`
- `plink` (1.9)
- `flashpca`
- `bcftools`
## Input data and toolpaths
The repository is equipped with a default mock-up input file `/data/pg1kg30x6.bisnv.nocompneighbor2.annot2.simple.1M.txt` which allows execution of steps 01-02 without any modifications.

> **NOTE** To run steps 03-07 substantial modification of Input data is required. Contact the authors for assistance in adapting the workflow to handle other input data (e.g. different features extracted from VCF files).

