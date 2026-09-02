> ## SUPERSEDED — DO NOT REUSE
>
> The analysis in this release was **withdrawn**. Errors were identified in dataset
> attribution, condition mapping, integer coercion, biological-unit assignment, and
> anatomical matching. Most consequentially, GSE116222 — the Parikh colonic
> *epithelial-isolation* experiment, not a multi-lineage atlas — was used as a
> multi-lineage deconvolution reference, and repeated biopsies from the same
> participant were analysed as independent observations.
>
> **The descriptions below are preserved exactly as they stood and contain those
> errors.** They are retained so that the corrections can be verified independently.
>
> The corrected analysis is release **`v2.0.0`**. A dated correction and deviation
> record accompanies the corrected manuscript as Additional file 2.

---
# SIRT6 — Epithelial-enriched inhibitor of the IBD inflammasome axis

Code, results, and figures for the manuscript:

> **SIRT6 is an epithelial-enriched inhibitor of the IBD inflammasome axis, with MAP3K7 (TAK1) as the most direct mRNA-level target: a single-cell, bulk, and external-validation analysis**

Pre-registered analysis plan (OSF): **DOI 10.17605/OSF.IO/AB5MY**

---

## Contents

```
github_upload/
├── README.md                       ← this file
├── scripts/                        28 .R / .Rmd analysis scripts
├── results/                        36 .csv result tables
├── figures/                        44 .png / .jpg figures (main + supplementary + QC)
└── preregistration/                3 PDFs (S1 + Amendments 01, 02)
```

| Folder | Files | Size | Description |
|--------|-------|------|-------------|
| `scripts/`        | 28 | 0.18 MB | All R / R Markdown analysis code |
| `results/`        | 36 | 1.61 MB | Per-stage CSV result tables |
| `figures/`        | 44 | 8.4 MB | PNG / JPG figures (Main Fig 1–5, Supplementary, IHC) |
| `preregistration/` | 3  | 0.55 MB | Pre-registered analysis plan + two amendments (PDF) |

---

## `scripts/` — Analysis code (28 files)

File-name prefix indicates the analytical stage:

| Prefix | Stage | Files |
|--------|-------|-------|
| `pilo_` | Stage 1 — GSE193677 bulk pilot | `01_pilot_GSE193677_v2.R`, package install scripts, `02_pilot_report.Rmd` |
| `scrn_` | Stage 2 — Single-cell + deconvolution | Smillie Seurat, MuSiC, HPA IHC, CIBERSORTx export, RELA targets, pseudobulk, Bisque |
| `exte_` | External validation + meta-analysis | GSE235236, GSE75214, 3-cohort meta-analysis, Fig 5D three-layer |
| `manu_` | Manuscript supplementary | metafor install, ppcor sensitivity, partial-Spearman triple validation, bootstrap CI |

Key scripts:
- `pilo_01_pilot_GSE193677_v2.R` — DESeq2 vst → SIRT6 Spearman → WGCNA on n = 2,483 cohort.
- `scrn_01_smillie_seurat_analysis.R` — Smillie 2019 scRNA SIRT6 cell-type localization.
- `scrn_02_music_deconvolution.R` — MuSiC cell fractions + partial Spearman.
- `scrn_06_pseudobulk_per_donor.R` — RELA paradox dropout-artifact control.
- `scrn_07_bisque_deconvolution.R` — Bisque cross-validation (replaces CIBERSORTx per S1-A01).
- `exte_01_GSE235236_validation.R` — RNA-Seq validation cohort #1.
- `exte_03_GSE75214_validation.R` — Microarray validation cohort #2 (added per S1-A02).
- `exte_04_3cohort_meta_analysis.R` — Random-effects REML pooled meta-analysis (k = 12 per target).
- `manu_supp_bootstrap_ci.R` — 1,000-resample bootstrap 95% CIs for partial Spearman.

Environment requirements: R 4.6.0+; Bioconductor 3.19+; Seurat 5+; metafor; BisqueRNA; ppcor; psych; clusterProfiler; SingleR + celldex (HumanPrimaryCellAtlas); WGCNA; DESeq2; org.Hs.eg.db; hugene10sttranscriptcluster.db (for GSE75214 microarray).

---

## `results/` — CSV tables (36 files)

| Prefix | Stage | Selected files |
|--------|-------|---------------|
| `pilo_` | GSE193677 pilot | `01_sirt6_correlation_table.csv`, `02_sirt6_module_hub_top30.csv`, `03_GO_BP_sirt6_module.csv` |
| `scrn_` | scRNA + deconvolution | `01_SIRT6_by_cell_type.csv`, `03_SIRT6_corr_within_celltype.csv`, `04_MuSiC_cell_fractions.csv`, `05_partial_correlation.csv`, `12_pseudobulk_per_donor_corr.csv`, `14_ppcor_sensitivity.csv`, `15_partial_triple_validation.csv`, `16_bootstrap_partial_CI.csv`, `17_Bisque_cell_fractions.csv`, `18_MuSiC_Bisque_agreement.csv` |
| `ext_GSE235236_` | Validation cohort #1 | `01_sirt6_correlation.csv`, `04_discovery_vs_validation.csv`, `05_meta_effect_sizes.csv`, `06_meta_analysis_pooled.csv`, `07_meta_3layer_pooled.csv`, `08_3cohort_effect_sizes.csv`, `10_3cohort_meta_pooled.csv`, `11_direction_consistency_per_target.csv` |
| `ext_GSE75214_` | Validation cohort #2 | `00_pdata.csv`, `01_sirt6_correlation.csv`, `02_sirt6_module_hub_top30.csv`, `03_disc_vs_GSE75214.csv` |

---

## `figures/` — Visual outputs (44 files)

Main figures referenced in the manuscript (15) plus supplementary, QC, and HPA IHC (29).

Selected key files:
- `Figure_1B_correlation_heatmap.png` — Main Fig 1B (bulk SIRT6 vs target heatmap).
- `Figure_3ABC_scRNA_SIRT6_localization.png` — Main Fig 3A–C (Smillie scRNA UMAP / dot plot).
- `Figure_3D_HPA_IHC_colon_donor1~3.jpg` — Main Fig 3D (Human Protein Atlas SIRT6 IHC).
- `Figure_4B_marginal_vs_partial.png` — Main Fig 4B (cell-composition-adjusted partial correlation).
- `Figure_5C_3cohort_meta_forest_84effects.png` — Main Fig 5C (random-effects forest plot, k = 12 per target across 3 cohorts).
- `Figure_5D_3cohort_pooled_sensitivity.png` — Main Fig 5D (3-cohort vs RNA-Seq-only sensitivity).
- `Supp_Figure_5A_MuSiC_vs_Bisque_scatter.png` — Bisque cross-validation per S1-A01.
- `Supp_Figure_SX_bootstrap_partial_CI.png` — 1,000-resample bootstrap 95% CIs.

For a complete index see the figure-naming convention used throughout: `Figure_N{A,B,C,...}_<description>` for main figures, `Supp_Figure_N_<description>` for supplementary, `Supp_QC_<description>` for QC-only.

---

## `preregistration/` — Pre-registered analysis plan (3 PDFs)

All deposited under OSF DOI **10.17605/OSF.IO/AB5MY** (https://osf.io/ab5my/).

| File | Frozen on | What it commits |
|------|-----------|-----------------|
| `S1_Pre-registered_Analysis_Plan.pdf` | 2026-05-21 | Stage 1 go/no-go rule, Stage 2 three-layer triangulation framework, multiple-testing handling (triangulation, no cross-layer correction; Lawlor et al. 2016), MuSiC↔CIBERSORTx Spearman r ≥ 0.7 agreement threshold. |
| `S1_Amendment_01_CIBERSORTx_to_Bisque.pdf` | 2026-05-21 21:49 KST | Methodological substitution: CIBERSORTx (Stanford server confirmed unavailable) → Bisque (R package BisqueRNA, Jew et al. 2020 Nat Commun). All thresholds and analytical commitments preserved verbatim. |
| `S1_Amendment_02_GSE75214_validation.pdf` | 2026-05-27 | Hypothesis-confirming addition of GSE75214 (n = 194 microarray, Vancamelbeke et al. 2017) as second external validation cohort to strengthen the meta-analytic pool (k = 8 → k = 12). All pre-registered statistical procedures, decision rules, and thresholds preserved. |

---

## Data sources (not included in this repository)

All raw and processed expression data are publicly available:

| Accession | Cohort | Platform | n |
|-----------|--------|----------|---|
| GSE193677 | Mt Sinai Crohn's and Colitis Registry (discovery) | Illumina HiSeq 2500 RNA-Seq | 2,483 mucosal biopsies |
| GSE116222 | Smillie et al. 2019 (single-cell) | 10× Chromium scRNA-Seq | 11,166 colonic cells |
| GSE235236 | External validation #1 | Illumina HiSeq 3000 RNA-Seq + RSEM | 56 sigmoid colon biopsies |
| GSE75214 | Vancamelbeke et al. 2017 (external validation #2) | Affymetrix Human Gene 1.0 ST microarray | 194 colonic/ileal biopsies |
| Human Protein Atlas | SIRT6 antibody HPA071776 | IHC images (colon, rectum, duodenum, pancreas, liver) | 26 images |

---

## Reproducibility

1. Install R 4.6.0+ and required packages (see comments at the top of each script; CRAN + Bioconductor + GitHub for MuSiC).
2. Download raw data from the GEO accessions listed above into the directory structure assumed by the scripts (paths are documented at the top of each `.R` file).
3. Run scripts in order: Stage 1 (`pilo_*`) → Stage 2 (`scrn_*`) → external validation (`exte_*`). Each script saves intermediate `.RDS` artifacts that downstream scripts read.
4. To re-run only the meta-analysis (which is the heaviest step depending on prior outputs), run `exte_04_3cohort_meta_analysis.R` after Stage 1 and the two external-validation scripts have completed.

---

## Citation

If you use any code or results from this repository, please cite:

- The main manuscript: [to be inserted once published]
- The pre-registered analysis plan: **OSF DOI 10.17605/OSF.IO/AB5MY**
- Underlying data: GSE193677 (Argmann et al.), GSE116222 (Smillie et al. 2019 Cell), GSE235236, GSE75214 (Vancamelbeke et al. 2017 Inflamm Bowel Dis)
- Methods: MuSiC (Wang et al. 2019 Nat Commun), Bisque (Jew et al. 2020 Nat Commun), Lawlor et al. 2016 Int J Epidemiol (triangulation framework), metafor (Viechtbauer 2010 J Stat Softw)

---

## Contact

Bang Changseok (방창석)
Project: SIRT6-IBD, 2026
For questions: [to be inserted]
