# ==============================================================================
# Stage 2 supplementary — Pseudobulk per donor 분석
# 목적: scRNA dropout artifact 차단으로 RELA paradox 진위 판정
# 방법: 각 donor의 broad cell type별 (특히 epithelial) raw count 합산 → pseudobulk
#       → SIRT6 vs RELA / 인플라마좀 target Spearman
# 데이터: Smillie GSE116222 (3 donors × 3 conditions = 9 sample / cell type)
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(Matrix)
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(ggpubr); library(pheatmap); library(scuttle); library(DESeq2)
})

set.seed(42)

# ==============================================================================
# 1. Load annotated Seurat
# ==============================================================================
cat(">>> [1/5] Load Smillie Seurat\n")
seu <- readRDS("rds/smillie_seurat_annotated.rds")
cat(sprintf("   %d cells × %d genes\n", ncol(seu), nrow(seu)))

# Define pseudobulk grouping: donor + condition + broad_type
# 9 donor-condition × 3-4 cell types = ~27-36 pseudobulk samples
seu$pb_group <- paste(seu$donor, seu$condition, seu$broad_type, sep = "_")
cat("   pseudobulk groups (donor × condition × broad_type):\n")
group_sizes <- table(seu$pb_group)
cat(sprintf("   total groups: %d, sizes range %d~%d\n",
            length(group_sizes), min(group_sizes), max(group_sizes)))

# Drop groups with < 10 cells (unreliable)
keep_groups <- names(group_sizes)[group_sizes >= 10]
cat(sprintf("   groups with ≥10 cells: %d\n", length(keep_groups)))
seu$pb_keep <- seu$pb_group %in% keep_groups

# ==============================================================================
# 2. Pseudobulk = sum of raw counts per group
# ==============================================================================
cat(">>> [2/5] Build pseudobulk matrix\n")
counts_raw <- GetAssayData(seu, assay = "RNA", layer = "counts")
keep_cells <- which(seu$pb_keep)

# aggregateAcrossCells from scuttle
sce <- SingleCellExperiment(
  assays = list(counts = counts_raw[, keep_cells]),
  colData = DataFrame(
    cell_id = colnames(seu)[keep_cells],
    donor = seu$donor[keep_cells],
    condition = seu$condition[keep_cells],
    broad_type = seu$broad_type[keep_cells],
    pb_group = seu$pb_group[keep_cells]
  )
)
pb <- aggregateAcrossCells(sce, ids = colData(sce)$pb_group, statistics = "sum")
cat(sprintf("   pseudobulk: %d genes × %d groups\n",
            nrow(pb), ncol(pb)))

# pb_meta
pb_meta <- as.data.frame(colData(pb)) %>%
  mutate(
    donor      = sub("_.*", "", pb_group),
    rest       = sub("^[^_]+_", "", pb_group),
    broad_type = sub(".*_", "", rest),
    condition  = sub("_[^_]+$", "", rest)
  ) %>% dplyr::select(pb_group, donor, condition, broad_type, ncells)
print(pb_meta)

# ==============================================================================
# 3. Normalize (DESeq2 vst, design ~broad_type+condition)
# ==============================================================================
cat(">>> [3/5] DESeq2 vst on pseudobulk\n")
dds_pb <- DESeqDataSetFromMatrix(
  countData = assay(pb, "counts"),
  colData = pb_meta,
  design = ~ broad_type + condition
)
dds_pb <- estimateSizeFactors(dds_pb)
vsd_pb <- vst(dds_pb, blind = FALSE)
expr_pb <- assay(vsd_pb)

# ==============================================================================
# 4. SIRT6 vs target Spearman — *separate per cell type*
# ==============================================================================
cat(">>> [4/5] Spearman SIRT6 vs targets (per broad_type)\n")
targets <- c("NLRP3", "CASP1", "IL1B", "PYCARD", "MAP3K7", "RELA", "FOXC1")
rela_targets <- c("IL6", "TNF", "CCL2", "CXCL8", "CXCL1", "ICAM1",
                  "NFKBIA", "TNFAIP3", "PTGS2", "BCL2", "BIRC3")
all_targets <- intersect(c(targets, rela_targets), rownames(expr_pb))

results <- list()
for (ct in unique(pb_meta$broad_type)) {
  cols_ct <- which(pb_meta$broad_type == ct)
  if (length(cols_ct) < 5) next
  s <- expr_pb["SIRT6", cols_ct]
  for (tg in all_targets) {
    t <- expr_pb[tg, cols_ct]
    if (sd(t) == 0) next
    ct_res <- cor.test(s, t, method = "spearman", exact = FALSE)
    results[[length(results) + 1]] <- data.frame(
      broad_type = ct, target = tg, n_samples = length(cols_ct),
      rho = unname(ct_res$estimate), p = ct_res$p.value
    )
  }
}
res_df <- bind_rows(results) %>%
  group_by(broad_type) %>%
  mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
  arrange(broad_type, p)

write.csv(res_df, "results/12_pseudobulk_per_donor_corr.csv", row.names = FALSE)
cat("Pseudobulk per-donor correlation:\n"); print(res_df, n = 50)

# ==============================================================================
# 5. Heatmap + 비교 (single-cell vs pseudobulk)
# ==============================================================================
cat(">>> [5/5] Visualization\n")

# Heatmap matrix
m <- res_df %>% dplyr::select(broad_type, target, rho) %>%
  pivot_wider(names_from = target, values_from = rho) %>%
  column_to_rownames("broad_type") %>% as.matrix()

# 7개 core target만 (Stage 1과 비교)
core <- intersect(c("CASP1", "NLRP3", "MAP3K7", "RELA", "PYCARD", "IL1B", "FOXC1"),
                  colnames(m))
png("figures/17_pseudobulk_core_targets.png", 1100, 350, res = 130)
print(pheatmap(m[, core, drop = FALSE],
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks = seq(-1, 1, length.out=101),
  main = "Pseudobulk per donor — SIRT6 vs inflammasome targets (Spearman ρ)",
  display_numbers = TRUE, fontsize_number = 9,
  cluster_rows = FALSE, cluster_cols = TRUE))
dev.off()

# RELA target panel
rela_in_m <- intersect(rela_targets, colnames(m))
if (length(rela_in_m) >= 3) {
  png("figures/18_pseudobulk_RELA_targets.png", 1200, 350, res = 130)
  print(pheatmap(m[, rela_in_m, drop = FALSE],
    color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    breaks = seq(-1, 1, length.out=101),
    main = "Pseudobulk per donor — SIRT6 vs RELA target genes",
    display_numbers = TRUE, fontsize_number = 9,
    cluster_rows = FALSE, cluster_cols = TRUE))
  dev.off()
}

# ==============================================================================
# 6. Comparison table: scRNA single-cell vs pseudobulk
# ==============================================================================
cat("\n>>> Single-cell vs Pseudobulk comparison\n")
sc_file <- "results/03_SIRT6_corr_within_celltype.csv"
if (file.exists(sc_file)) {
  sc_df <- read.csv(sc_file) %>%
    dplyr::select(broad, target, rho_sc = rho)
  cmp <- res_df %>%
    dplyr::select(broad_type, target, rho_pb = rho) %>%
    inner_join(sc_df, by = c("broad_type" = "broad", "target")) %>%
    mutate(consistent = sign(rho_sc) == sign(rho_pb))
  print(cmp)
  write.csv(cmp, "results/13_scRNA_vs_pseudobulk_comparison.csv", row.names = FALSE)
  cat(sprintf("\n   Direction consistent: %d / %d\n",
              sum(cmp$consistent, na.rm = TRUE), nrow(cmp)))
}

cat("\n=== Pseudobulk DONE ===\n")
cat("  results/12_pseudobulk_per_donor_corr.csv\n")
cat("  results/13_scRNA_vs_pseudobulk_comparison.csv\n")
cat("  figures/17_pseudobulk_core_targets.png\n")
cat("  figures/18_pseudobulk_RELA_targets.png\n")
