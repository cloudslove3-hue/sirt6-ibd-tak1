# ==============================================================================
# Smillie 2019 scRNA — SIRT6 cell-type expression 분석
# 데이터: GSE116222 (11,174 cells, 3 donors × {Healthy, UC inflamed, UC non-inflamed})
#         원 paper의 validation cohort. Discovery cohort(SCP259, 366K cells)는
#         별도 등록 필요. 본 분석은 validation cohort로 1차 검증.
# 출력: SIRT6 cell-type localization → 1차 paper Fig 3 prototype
# 환경: R 4.6.0, Seurat 5+, C:/Rlibs, C:/Temp/R
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))

setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")
for (d in c("results", "figures", "rds")) dir.create(d, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(viridis); library(SingleR); library(celldex)
  library(SingleCellExperiment); library(scuttle)
  library(data.table); library(stringr); library(tibble)
})

set.seed(42)

# ==============================================================================
# 1. Load expression matrix
# ==============================================================================
cat(">>> [1/8] Loading GSE116222 expression matrix\n")
expr_file <- "data/GSE116222_Expression_matrix.txt.gz"
# 파일 구조: header=11,175 cells, 데이터행은 gene_id + 11,175 values (총 11,176 토큰)
# GSE193677과 동일하게 explicit header/body 분리
hdr <- readLines(gzfile(expr_file), n = 1)
sample_ids <- strsplit(hdr, "\t")[[1]]
sample_ids <- gsub('^"|"$', '', sample_ids)
cat(sprintf("   header cells: %d, 첫 3개: %s\n",
            length(sample_ids), paste(head(sample_ids, 3), collapse = " | ")))

expr_dt <- fread(expr_file, sep = "\t", header = FALSE, skip = 1,
                 data.table = FALSE, quote = "")
cat(sprintf("   body dim: %d x %d\n", nrow(expr_dt), ncol(expr_dt)))
stopifnot(ncol(expr_dt) == length(sample_ids) + 1)
gene_ids <- expr_dt[, 1]
expr_mat <- as.matrix(expr_dt[, -1])
storage.mode(expr_mat) <- "integer"
rownames(expr_mat) <- gene_ids
colnames(expr_mat) <- sample_ids
cat(sprintf("   matrix: %d genes x %d cells. 첫 cell: %s\n",
            nrow(expr_mat), ncol(expr_mat), head(sample_ids, 1)))
rm(expr_dt); gc()

# Sample annotation from cell barcode suffix
# A/B/C = donor 1/2/3, 1/2/3 = Healthy / UC inflamed / UC non-inflamed
barcodes <- colnames(expr_mat)
suffix <- sub(".*-", "", barcodes)
donor <- substr(suffix, 1, 1)
cond_code <- substr(suffix, 2, 2)
condition <- factor(
  ifelse(cond_code == "1", "Healthy",
   ifelse(cond_code == "2", "UC_inflamed", "UC_noninflamed")),
  levels = c("Healthy", "UC_noninflamed", "UC_inflamed")
)
cat("   Sample distribution:\n")
print(table(donor, condition))

# Sparse matrix
expr_sparse <- as(expr_mat, "CsparseMatrix")
rm(expr_mat, expr); gc()

# ==============================================================================
# 2. Build Seurat object + QC
# ==============================================================================
cat(">>> [2/8] Build Seurat + QC\n")
seu <- CreateSeuratObject(counts = expr_sparse, project = "Smillie2019_validation",
                          min.cells = 3, min.features = 200)
seu$donor <- donor[match(colnames(seu), barcodes)]
seu$condition <- condition[match(colnames(seu), barcodes)]
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

p_qc <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0)
ggsave("figures/01_QC_violin.png", p_qc, width = 12, height = 4, dpi = 120)

# Filter: cells with reasonable QC
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                            percent.mt < 25)
cat(sprintf("   QC pass: %d cells\n", ncol(seu)))

# ==============================================================================
# 3. Normalize + PCA + UMAP + clustering (캐싱)
# ==============================================================================
cache_clust <- "rds/seu_clustered.rds"
if (file.exists(cache_clust)) {
  cat(">>> [3/8] Loading cached clustered Seurat\n")
  seu <- readRDS(cache_clust)
} else {
  cat(">>> [3/8] Normalize + dimensionality reduction\n")
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2500, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

  p_elbow <- ElbowPlot(seu, ndims = 50)
  ggsave("figures/02_PCA_elbow.png", p_elbow, width = 7, height = 4, dpi = 120)

  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  saveRDS(seu, cache_clust)
}

# ==============================================================================
# 4. SingleR cell type annotation
# ==============================================================================
cat(">>> [4/8] SingleR cell-type annotation\n")
cache_sr <- "rds/seu_singler.rds"
if (file.exists(cache_sr)) {
  cat("   Loading cached SingleR\n")
  seu <- readRDS(cache_sr)
} else {
  ref <- celldex::HumanPrimaryCellAtlasData()
  sce <- as.SingleCellExperiment(seu)
  # de.method = "classic" works without scrapper dependency
  sr <- SingleR(test = sce, ref = ref, labels = ref$label.main,
                de.method = "classic")
  seu$singler_main <- sr$labels
  sr_fine <- SingleR(test = sce, ref = ref, labels = ref$label.fine,
                     de.method = "classic")
  seu$singler_fine <- sr_fine$labels
  saveRDS(seu, cache_sr)
}

# Broad category mapping (epithelial / immune / stromal)
broad_map <- function(lab) {
  case_when(
    grepl("Epithelial|Enterocyte|Goblet|Paneth|Tuft|Keratinocyte|Hepatocyte",
          lab, ignore.case = TRUE) ~ "Epithelial",
    grepl("T_cell|T cell|B_cell|B cell|NK|Macrophage|Monocyte|DC|Dendritic|Neutrophil|Mast|Plasma",
          lab, ignore.case = TRUE) ~ "Immune",
    grepl("Fibroblast|Endothelial|Smooth muscle|Stromal|Pericyte",
          lab, ignore.case = TRUE) ~ "Stromal",
    TRUE ~ "Other"
  )
}
seu$broad_type <- broad_map(seu$singler_main)
cat("   Broad type distribution:\n"); print(table(seu$broad_type))
cat("   Main label distribution:\n"); print(sort(table(seu$singler_main), decreasing = TRUE))

# ==============================================================================
# 5. SIRT6 expression by cell type
# ==============================================================================
cat(">>> [5/8] SIRT6 expression analysis\n")
key_genes <- c("SIRT6", "NLRP3", "CASP1", "IL1B", "PYCARD",
               "MAP3K7", "RELA", "FOXC1",
               "HNF4A", "EPCAM", "VIL1", "MUC2", "LYZ",  # epithelial markers
               "CD3D", "CD8A", "CD4", "FCGR3A", "CD68",   # immune markers
               "COL1A1", "VWF", "ACTA2")                  # stromal markers
key_genes <- intersect(key_genes, rownames(seu))

# Summary stats per cell type
expr_norm <- GetAssayData(seu, assay = "RNA", layer = "data")
sirt6_vals <- expr_norm["SIRT6", ]

sirt6_by_type <- data.frame(
  cell = colnames(seu),
  SIRT6 = sirt6_vals,
  broad = seu$broad_type,
  main = seu$singler_main,
  cluster = seu$seurat_clusters,
  condition = seu$condition,
  donor = seu$donor
)

# Cell type ranking by SIRT6 mean expression
sirt6_rank <- sirt6_by_type %>%
  group_by(main) %>%
  summarise(mean_SIRT6 = mean(SIRT6),
            median_SIRT6 = median(SIRT6),
            pct_expressing = mean(SIRT6 > 0) * 100,
            n_cells = n(), .groups = "drop") %>%
  arrange(desc(mean_SIRT6))

write.csv(sirt6_rank, "results/01_SIRT6_by_cell_type.csv", row.names = FALSE)
cat("   SIRT6 ranking by main cell type:\n"); print(sirt6_rank)

# Broad category SIRT6 stats by condition
sirt6_broad <- sirt6_by_type %>%
  group_by(broad, condition) %>%
  summarise(mean_SIRT6 = mean(SIRT6),
            pct_expressing = mean(SIRT6 > 0) * 100,
            n_cells = n(), .groups = "drop")
write.csv(sirt6_broad, "results/02_SIRT6_broad_by_condition.csv", row.names = FALSE)

# Save Seurat object
saveRDS(seu, "rds/smillie_seurat_annotated.rds")
saveRDS(sirt6_by_type, "rds/sirt6_per_cell.rds")

# ==============================================================================
# 6. UMAP + Dot plot + Violin (Fig 3 prototype components)
# ==============================================================================
cat(">>> [6/8] Fig 3 prototype components\n")

# 6-1. UMAP — broad cell type
p_umap_broad <- DimPlot(seu, group.by = "broad_type",
                        label = TRUE, repel = TRUE) +
  ggtitle("Broad cell type (SingleR + HPCA)") +
  theme(legend.position = "right")
ggsave("figures/03_UMAP_broad.png", p_umap_broad,
       width = 8, height = 6, dpi = 130)

# 6-2. UMAP — fine cell type
p_umap_main <- DimPlot(seu, group.by = "singler_main",
                       label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("Main cell type") +
  theme(legend.position = "none")
ggsave("figures/04_UMAP_main.png", p_umap_main,
       width = 9, height = 7, dpi = 130)

# 6-3. UMAP — SIRT6 expression
p_umap_sirt6 <- FeaturePlot(seu, features = "SIRT6",
                             cols = c("lightgrey", "darkred"),
                             pt.size = 0.4, order = TRUE) +
  ggtitle("SIRT6 expression") +
  theme(plot.title = element_text(face = "italic"))
ggsave("figures/05_UMAP_SIRT6.png", p_umap_sirt6,
       width = 7, height = 6, dpi = 130)

# 6-4. UMAP — condition
p_umap_cond <- DimPlot(seu, group.by = "condition",
                       cols = c("steelblue", "orange", "firebrick")) +
  ggtitle("Condition")
ggsave("figures/06_UMAP_condition.png", p_umap_cond,
       width = 8, height = 6, dpi = 130)

# 6-5. Dot plot — cell type × key genes
p_dot <- DotPlot(seu, features = key_genes,
                  group.by = "singler_main", scale = TRUE) +
  RotatedAxis() +
  scale_color_viridis_c(option = "magma") +
  ggtitle("Cell type × key gene expression (SIRT6 vs inflammasome vs markers)")
ggsave("figures/07_DotPlot_cell_x_gene.png", p_dot,
       width = 14, height = 7, dpi = 130)

# 6-6. Violin — SIRT6 by broad type, split by condition
p_vln_sirt6 <- VlnPlot(seu, features = "SIRT6", group.by = "broad_type",
                       split.by = "condition", pt.size = 0) +
  ggtitle("SIRT6 by broad cell type × condition")
ggsave("figures/08_Violin_SIRT6_broad.png", p_vln_sirt6,
       width = 9, height = 5, dpi = 130)

# 6-7. Combined Fig 3 prototype
fig3_prototype <- (p_umap_broad + p_umap_sirt6) /
                  (p_dot)
ggsave("figures/09_Fig3_prototype.png", fig3_prototype,
       width = 16, height = 11, dpi = 130)

# ==============================================================================
# 7. SIRT6 vs inflammasome co-expression (cell-level)
# ==============================================================================
cat(">>> [7/8] Cell-level SIRT6 vs target gene correlation (per cell type)\n")

targets <- c("NLRP3", "CASP1", "IL1B", "PYCARD", "MAP3K7", "RELA", "FOXC1")
targets <- intersect(targets, rownames(expr_norm))

cor_by_type <- list()
for (ct in unique(seu$broad_type)) {
  cells_ct <- which(seu$broad_type == ct)
  if (length(cells_ct) < 50) next
  s <- expr_norm["SIRT6", cells_ct]
  for (tg in targets) {
    t <- expr_norm[tg, cells_ct]
    # Drop double-zero pairs for cleaner signal
    nz <- (s > 0) | (t > 0)
    if (sum(nz) < 30) next
    ct_res <- cor.test(s[nz], t[nz], method = "spearman", exact = FALSE)
    cor_by_type[[length(cor_by_type) + 1]] <- data.frame(
      broad = ct, target = tg, n_cells_nz = sum(nz),
      rho = unname(ct_res$estimate), p = ct_res$p.value
    )
  }
}
cor_celltype <- bind_rows(cor_by_type) %>%
  group_by(broad) %>%
  mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
  arrange(broad, p)
write.csv(cor_celltype, "results/03_SIRT6_corr_within_celltype.csv",
          row.names = FALSE)
cat("   SIRT6 vs targets correlation within cell types (cell-level):\n")
print(cor_celltype %>% arrange(target, broad))

# Visualize
cor_wide <- cor_celltype %>%
  select(broad, target, rho) %>%
  pivot_wider(names_from = target, values_from = rho) %>%
  column_to_rownames("broad") %>% as.matrix()
p_celltype_cor <- pheatmap::pheatmap(cor_wide,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-0.4, 0.4, length.out = 101),
  main = "SIRT6 vs inflammasome — within cell type (Spearman ρ)",
  display_numbers = TRUE, fontsize_number = 9,
  filename = "figures/10_corr_within_celltype.png",
  width = 8, height = 4)

# ==============================================================================
# 8. Pseudobulk for MuSiC reference (다음 script C에서 사용)
# ==============================================================================
cat(">>> [8/8] Save reference matrix for MuSiC deconvolution\n")
sce_for_music <- as.SingleCellExperiment(seu)
colData(sce_for_music)$cellType <- seu$broad_type  # broad categories
colData(sce_for_music)$cellType_fine <- seu$singler_main
colData(sce_for_music)$subjectID <- paste0(seu$donor, "_", seu$condition)

saveRDS(sce_for_music, "rds/smillie_sce_for_music.rds")
cat("   SCE saved for MuSiC reference: rds/smillie_sce_for_music.rds\n")

cat("\n=== Smillie scRNA analysis DONE ===\n")
cat("핵심 출력:\n")
cat("  - results/01_SIRT6_by_cell_type.csv (cell type별 SIRT6 평균 발현 순위)\n")
cat("  - results/02_SIRT6_broad_by_condition.csv\n")
cat("  - results/03_SIRT6_corr_within_celltype.csv (cell type 내 co-expression)\n")
cat("  - figures/03~10.png (UMAP, dot plot, violin, heatmap)\n")
cat("  - figures/09_Fig3_prototype.png ★ 1차 paper Fig 3 prototype\n")
cat("  - rds/smillie_seurat_annotated.rds (Seurat 객체, 후속 분석용)\n")
cat("  - rds/smillie_sce_for_music.rds (MuSiC deconvolution용)\n")
