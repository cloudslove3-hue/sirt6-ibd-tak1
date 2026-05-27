# Figure files — submission package mapping

총 39개 figure (15 main + 24 supplementary/QC) — manuscript v3 인용 순서대로 rename.

## Main Figures (15 files)

| Manuscript | File | 설명 |
|------------|------|------|
| Fig 1B | `Figure_1B_correlation_heatmap.png` | SIRT6 vs 7 target Spearman heatmap (8 subgroup × 7 target), n=2,483 |
| Fig 1C | `Figure_1C_UC_inflamed_scatter.png` | UC_inflamed scatter (n=293, 7 panel) |
| Fig 2A | `Figure_2A_WGCNA_softthresh.png` | WGCNA soft-thresholding (power 14 선택) + connectivity |
| Fig 2C | `Figure_2C_GO_BP_top10.png` | SIRT6 module GO BP top 10 (epithelial barrier) |
| Fig 3A-C | `Figure_3ABC_scRNA_SIRT6_localization.png` | Smillie scRNA: UMAP broad + SIRT6 expression + dotplot |
| Fig 3D | `Figure_3D_HPA_IHC_colon_donor1.jpg` | HPA IHC SIRT6 colon, donor 1 (crypt + surface enterocyte) |
| Fig 3D | `Figure_3D_HPA_IHC_colon_donor2.jpg` | HPA IHC SIRT6 colon, donor 2 |
| Fig 3D | `Figure_3D_HPA_IHC_colon_donor3.jpg` | HPA IHC SIRT6 colon, donor 3 |
| Fig 3E | `Figure_3E_within_celltype_heatmap.png` | Within-cell-type Spearman heatmap (Epi/Imm × 7 target) |
| Fig 4A | `Figure_4A_MuSiC_cell_fractions.png` | MuSiC cell fraction by disease × inflammation |
| Fig 4B | `Figure_4B_marginal_vs_partial.png` | ★ Marginal vs partial correlation bar plot (4 subgroup × 7 target) |
| Fig 5A | `Figure_5A_GSE235236_heatmap.png` | GSE235236 validation heatmap |
| Fig 5B | `Figure_5B_discovery_vs_validation_scatter.png` | Discovery vs Validation ρ scatter (28쌍) |
| Fig 5C | `Figure_5C_meta_forest_56effects.png` | ★ Random-effects forest plot (56 effect sizes) |
| Fig 5D | `Figure_5D_meta_pooled_3layer.png` | ★ Per-target pooled ρ — Discovery/Validation/All-studies 3 layer |

## Supplementary Figures (24 files)

### Documented in manuscript v3
| Supp | File | 설명 |
|------|------|------|
| S Fig 1 (rectum) | `Supp_Figure_1_HPA_IHC_rectum.jpg` | HPA IHC SIRT6 rectum |
| S Fig 1 (duodenum) | `Supp_Figure_1_HPA_IHC_duodenum.jpg` | HPA IHC SIRT6 duodenum |
| S Fig 1 (pancreas) | `Supp_Figure_1_HPA_IHC_pancreas.jpg` | HPA IHC SIRT6 pancreas |
| S Fig 1 (liver) | `Supp_Figure_1_HPA_IHC_liver.jpg` | HPA IHC SIRT6 liver |
| S Fig 2A | `Supp_Figure_2A_pseudobulk_core_targets.png` | Pseudobulk per donor — 7 core target heatmap |
| S Fig 2B | `Supp_Figure_2B_pseudobulk_RELA_targets.png` | Pseudobulk per donor — RELA target heatmap |
| S Fig 3A | `Supp_Figure_3A_RELA_targets_marginal.png` | RELA target gene panel — marginal correlation |
| S Fig 3B | `Supp_Figure_3B_RELA_targets_partial.png` | RELA target gene panel — partial correlation |
| S Fig 5A | `Supp_Figure_5A_MuSiC_vs_Bisque_scatter.png` | ★ Bisque cross-validation scatter (Epi + Imm) |
| S Fig 5B | `Supp_Figure_5B_BlandAltman_MuSiC_Bisque.png` | Bland-Altman MuSiC vs Bisque |
| S Fig SX | `Supp_Figure_SX_bootstrap_partial_CI.png` | ★ Bootstrap 95% CI (28 combinations) |

### QC / supporting (manuscript에 직접 인용 안 됨, repo 보존용)
- `Supp_QC_GSE193677_library_size.png` — GSE193677 library size distribution
- `Supp_KEGG_top10.png` — SIRT6 module KEGG (q > 0.3, supplementary)
- `Supp_QC_Smillie_violin.png` — Smillie scRNA QC violin
- `Supp_QC_Smillie_PCA_elbow.png` — Smillie PCA elbow
- `Supp_QC_Smillie_UMAP_broad.png` — Smillie UMAP broad cell type
- `Supp_QC_Smillie_UMAP_main.png` — Smillie UMAP main cell type
- `Supp_QC_Smillie_UMAP_SIRT6.png` — Smillie UMAP SIRT6 expression
- `Supp_QC_Smillie_UMAP_condition.png` — Smillie UMAP by condition
- `Supp_QC_Smillie_dotplot.png` — Smillie dot plot
- `Supp_QC_Smillie_violin_SIRT6.png` — SIRT6 violin by broad type
- `Supp_QC_partial_heatmap_marginal.png` — Partial heatmap (marginal)
- `Supp_QC_partial_heatmap_partial.png` — Partial heatmap (epi-adjusted)
- `Supp_QC_ppcor_algebraic_equivalence.png` — ppcor algebraic equivalence scatter (Pearson r=1.0000)

## 원본 위치 (재생성용)

| 원본 폴더 | 파일 수 |
|---------|--------|
| `pilot_analysis/figures/` | 6개 |
| `scrna_analysis/figures/` | 20+개 |
| `scrna_analysis/hpa/images/` | 26개 |
| `external_validation/figures/` | 5개 |

## CMGH 제출 시 figure 처리

CMGH submission portal은 보통 figure를 별도 high-res 파일 (TIFF, EPS, 또는 PDF)로 업로드.
현재는 PNG/JPG (~150-300 DPI). 제출 직전:
1. ggplot 객체 다시 호출 → `ggsave(..., dpi = 300, device = "tiff")`로 재렌더 권장
2. 모든 Main Figure를 300 DPI TIFF로 변환
3. Supplementary는 PNG 그대로 허용
