#  Mitigating batch effects in multicenter whole genome sequencing with duplicate samples and machine learning

[![DOI](https://zenodo.org/badge/DOI/ZENODO_DOI_HERE.svg)](https://doi.org/ZENODO_DOI_HERE)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code used in **Mitigating batch effects in multicenter whole genome sequencing with duplicate samples and machine learning** ([journal link/DOI once available]). The paper presents a proof-of-concept approach of WGS data QC utilizing duplicate samples and machine learning. This repository is intended to document and accompany the methods described in the paper, and is **not** a bioinformatics tool itself.

Due to the use of both publicly accessible and private WGS datasets, full end-to-end reproducibility is not possible. See the [Data Availability](#data-availability) section for details.

---


## Repository Structure

The repository is organised into three self-contained workflows, each corresponding to a QC strategy described in the paper. Within each workflow, scripts are numbered and intended to be executed sequentially. A description of the directory layout is provided below in the [Workflows](#workflows) section.

```
DuplicateSampleBatchEffectFilter/
├── README.md
├── requirements.txt
├── LociBased/
│   ├── README.md
│   ├── data/
│   ├── models/
│   ├── outputs/
│   └── scripts/
├── SampleBased/
│   ├── data/
│   ├── models/
│   ├── outputs/
│   └── scripts/
└── OnePair/
  ├── data/
  ├── models/
  ├── outputs/
  └── scripts/
```

---

## Requirements
- bcftools 1.18
- PLINK v1.90b6.21 64-bit (19 Oct 2020)
- PLINK v2.00a4.4LM AVX2 Intel (21 Jun 2023) 
- flashpca 2.0
- git lfs
- __Python__ 3.13.7
    - pandas 2.3.2
    - numpy 2.3.2
    - scipy 1.16.1
    - xgboost 3.0.4
    - scikit-learn 1.7.1
    - matplotlib 3.10.6

- __R__ 4.5.2
    - ggplot2 3.5.2
    - ggrepel 0.9.6
    - dplyr 1.1.4
    - CMplot 4.51.1

---

## Workflows

All three self-contained workflows follow the same sequential structure: 

### 1.  LociBased strategy — initial design of WGS data batch effect correction using multisample data ###
| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_orchestration.py` | run entire workflow step by step |
| 1 | `01_prepare_training_data.py` | preprocessing of data exported from VCF|
| 2 | `02_train_xgb_classifier.py.py` | XGBoost model from paired sampels |
| 3 | `03_predict_variants.py` | predict reliable/likely erroneous variants |
| 4 | `04_filter_bcf_file.py` | remove predicted varaints from multisample .bcf file |
| 5 | `05_bcf2pca.sh` |  prepare PCA from filtered .bcf file |
| 6 | `06_vec2pca.R` | plot PCA |
| 7 | `07_mhqq.R`| mock-up GWAS with Manhattan and QQ plots |

Workflow demonstrates initial strategy of paired sample utilization. Relative simplicity of this strategy (multisample data organozation) introduces noisy labels, which compromise overall accuracy. Utilizes 6 duplicate samples across two sequencing center data.

> **Note:** Steps 01-02 of this workflow can be run using publicly available data (small example data enabled by default). See [Data Availability](#data-availability) for download instructions. Steps 03-07 requires access to private WGS data that cannot be shared. Scripts are provided for transparency and methodological documentation only.

### 2. SampleBased strategy — batch effect correction with sample-by-sample processing and weak-supevision labeling
| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_orchestration.py` | run entire workflow step by step |
| 1 | `01_weak_sup_labeling.py` | labeling of data exported from VCF|
| 2 | `02_prepare_training_data.py` | prepare labeled data for ML |
| 3 | `03_train_xgb_classifier.py.py` | XGBoost model from paired sampel data |
| 4 | `04_predict_variants.py` | predict reliable/likely erroneous variants |
| 5 | `05_filter_bcf_files.py` | filter single-sample .bcf files one by one |
| 6 | `06_merge_bcf_files.py` | merge single-sampe .bcf files into one |
| 7 | `07_bcf2pca.sh` |  prepare PCA from filtered .bcf file |
| 8 | `08_vec2pca.R` | plot PCA |
| 9 | `09_mhqq.R`| mock-up GWAS with Manhattan and QQ plots |

Workflow adresses the noisy label issue of the loci-based strategy but increases complexity and processing time, improves filtering accuracy by utilizing 6 duplicate samples across two sequencing center data.

> **Note:** Steps 01-03 of this workflow can be run using publicly available data (small exampkle data enabled by default, no workflow modifications required). See [Data Availability](#data-availability) for download instructions.  Steps 04-09 requires access to private WGS data that cannot be shared. Scripts are provided for transparency and methodological documentation only.

### 3. OnePair strategy — one paired (duplicate) sample is enough for Batch effect filtering 

Essentialy a copy of Sample-based strategy with only one duplicate sample as training data. Demonstrates similar quality of batch effect diltering while requiring smaller training data.
The workflow has the same structure, input data, and script order as Sample-based strategy. 


---

## Running the Code
Clone repository and get the input data.
```bash
git clone https://github.com/Le-Crab/DuplicateSampleBatchEffectFilter
cd DuplicateSampleBatchEffectFilter
git lfs install
git lfs pull # get example subset of input data to run the workflows
```
Execute scripts within each workflow in numerical order. For example, for Loci-based stratgy:

```bash
# execuion of steps one by one
cd LociBased/scripts
python 01_prepare_training_data.py
python 02_train_xgb_classifier.py 

# run full workflow 
cd LociBased/scripts
python 00_orchestration.py
```
> **Note:**
Workflow steps after training the ML model require private data and cannot be fully executed without it. The scripts are included to fully document the methods described in the paper. For assistance in adapting the workflow to fit custom data, please contact the authors.

---

## Data Availability

### Public Data

| Dataset  |Link |
|---------|-------------------|
| Platinum genomes | https://www.ebi.ac.uk/ena/browser/view/PRJEB3381|
| 1000 Genomes Project| https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x |

### Private Data
| Dataset  |Link |
|---------|-------------------|
| Biobank Japan  | https://humandbs.dbcls.jp/en/hum0014-v32#JGAS000114 |
| Yaponesia |  available upon request to the paper correspondence author|

For full details, refer to the **Data Availability** statement in the paper.

---

## Citation

If you use this code, please cite:

```bibtex
@article{[citekey],
  author  = {[Authors]},
  title   = {[Paper Title]},
  journal = {Bioinformatics},
  year    = {[year]},
  doi     = {[DOI]}
}
```

This repository is archived at Zenodo: [https://doi.org/ZENODO_DOI_HERE](https://doi.org/ZENODO_DOI_HERE)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

Uladzislau Korzun — crab.cnx@gmail.com

For questions about the code, please open an [issue](../../issues) in this repository.