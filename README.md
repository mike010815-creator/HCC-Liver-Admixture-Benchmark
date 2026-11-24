# HCC Non-Tumor Liver RNA Admixture Benchmark & Single-Cell-Guided Gene Purification

[![R](https://img.shields.io/badge/R-4.6.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Data: Public](https://img.shields.io/badge/Data-TCGA%20%7C%20ICGC%20%7C%20GEO-blue)](#public-data-sources)

Official code repository and curated gene-set resources accompanying the manuscript:

> **Non-Tumor Liver RNA Admixture Systematically Distorts Immune Signatures in Hepatocellular Carcinoma: Multi-Modal Dissection and Single-Cell-Guided Gene Purification**  
> Shifeng Xie<shifengxie@stu.2026.jun.edu.cn>, Chengzhi Li<lichengzhi@jnu.edu.cn>

---

## Overview

Bulk transcriptomic immune and tumor microenvironment (TME) signatures in hepatocellular carcinoma (HCC) are frequently scored on surgical resections or core-needle biopsies containing variable proportions of inflamed adjacent liver parenchyma. This repository provides:

1. **Multi-modal benchmarking pipeline** across **439 matched tumor–adjacent pairs (878 profiles)** from TCGA-LIHC, ICGC-LIRI-JP (HCCDB18), and GSE14520, **110 healthy GTEx livers** (Toil-processed), **paired scRNA-seq** (GSE149614; 7 paired donors), **all 7 Visium spatial sections** (GSE238264), and **789 primary HCC survival records** (763 complete cases for multivariable Cox adjustment).
2. **Ready-to-use gene-set resources**:
   - **NHAS (22 genes)**: Zero-overlap hepatocyte-associated expression score derived strictly from GSE149614 (`resources/Supplementary_Data_S3_NHAS_22_Genes.csv`).
   - **Original 24 gene sets** (`resources/Supplementary_Data_S1_Original_24_GeneSets.gmt`) and **18 single-cell-purified (S3) gene sets** (`resources/Supplementary_Data_S2_S3_Purified_GeneSets.gmt`) with zero bulk data leakage.
3. **Statistical validation suite**: Patient-paired bootstrap 95% CIs ($B = 1{,}000$), zero-leakage uniform and expression-bin-matched random gene-deletion null controls ($R = 50$), cross-cohort external-reference subtraction (S2), non-monotonic NHAS quartile diagnostics, unified all-section spatial scoring, and multivariable Cox meta-analysis.

---

## Repository Structure

```text
.
├── README.md
├── scripts/
│   ├── 01_data_prep_*.R                         # Bulk, Toil GTEx, scRNA-seq, and Visium preprocessing
│   ├── 02_bulk_admixture_sim.R                  # Linear-scale 0%–50% admixture across ssGSEA / GSVA / Z-score
│   ├── 03_scrna_mechanism_and_nhas.R            # Dissociation-reweighted pseudo-bulk dissection & 22-gene NHAS
│   ├── 04_clinical_survival_audit.R             # 789-patient real-tumor NHAS co-variation & univariable Meta-Cox
│   ├── 05_correction_benchmark.R                # Head-to-head S1 / S2 / S3 deconfounding benchmark
│   ├── 07a_bootstrap_null_and_external_S2.R     # Patient bootstrap (B=1,000), matched null (R=50), external S2
│   ├── 07b_interaction_clinical_spatial.R       # NHAS quartile diagnostics, multivariable Cox (n=763), all-7 Visium
│   ├── 08b_polish_figs_4_to_6.R                 # Publication-grade Figures 4–6 (ggplot2 / patchwork)
│   └── 09_build_submission_package.R            # Read-only packaging & MD5 manifest generator
├── resources/
│   ├── Supplementary_Data_S1_Original_24_GeneSets.gmt
│   ├── Supplementary_Data_S2_S3_Purified_GeneSets.gmt
│   ├── Supplementary_Data_S3_NHAS_22_Genes.csv
│   └── Supplementary_Tables.xlsx                # Complete Tables S1–S15 (coverage, atlas, bootstrap, null, Cox)
└── env/
    ├── R_package_versions.csv
    └── R_sessionInfo.txt
