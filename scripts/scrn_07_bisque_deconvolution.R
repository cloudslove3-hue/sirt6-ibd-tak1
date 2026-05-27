# ==============================================================================
# Bisque reference-based deconvolution — CIBERSORTx 대체
# 사유: CIBERSORTx Stanford server confirmed down 2026-05-21 21:49 KST
# Method: BisqueRNA::ReferenceBasedDecomposition (Jew et al. 2020 Nat Commun)
# Reference: 동일한 Smillie GSE116222 broad cell type
# Mixture: 동일한 GSE193677 bulk (n=2,483)
# 출력: Bisque epithelial fraction → MuSiC과 cross-validation
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- 1. Bisque 설치 ----
# CRAN 첫 시도; 실패 시 GitHub
if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
  cat("Trying CRAN BisqueRNA...\n")
  tryCatch(
    install.packages("BisqueRNA", lib = "C:/Rlibs"),
    error = function(e) cat("CRAN failed:", conditionMessage(e), "\n")
  )
}
if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
  cat("Trying GitHub cozygene/bisque...\n")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", lib = "C:/Rlibs")
  tryCatch(
    remotes::install_github("cozygene/bisque", upgrade = "never",
                            lib = "C:/Rlibs", dependencies = TRUE),
    error = function(e) cat("GitHub failed:", conditionMessage(e), "\n")
  )
}
if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
  stop("BisqueRNA install failed via both CRAN and GitHub")
}

suppressPackageStartupMessages({
  library(BisqueRNA); library(Biobase); library(SingleCellExperiment)
  library(Matrix); library(dplyr); library(tidyr); library(tibble)
  library(data.table); library(org.Hs.eg.db); library(AnnotationDbi)
  library(matrixStats); library(ggplot2)
})
set.seed(42)

# ---- 2. Load Smillie reference (SCE) ----
cat(">>> [1/6] Load Smillie SCE\n")
sce <- readRDS("rds/smillie_sce_for_music.rds")
cat(sprintf("   Smillie: %d cells × %d genes\n", ncol(sce), nrow(sce)))

# Bisque는 ExpressionSet 형식 요구
sc_counts <- assay(sce, "counts")
sc_meta <- as.data.frame(colData(sce))
# Bisque 표준 컬럼명: cellType (이미 있음) + SubjectName
sc_meta$SubjectName <- sc_meta$subjectID

sc_eset <- ExpressionSet(
  assayData = as.matrix(sc_counts),
  phenoData = AnnotatedDataFrame(sc_meta[, c("cellType", "SubjectName")])
)
cat(sprintf("   SC ExpressionSet: %d genes × %d cells\n",
            nrow(sc_eset), ncol(sc_eset)))
cat("   cellType:\n"); print(table(sc_eset$cellType))

# ---- 3. Load GSE193677 bulk (raw counts) ----
cat(">>> [2/6] Load GSE193677 raw counts\n")
counts_file <- "../pilot_analysis/data/GSE193677_counts.txt.gz"
hdr <- readLines(gzfile(counts_file), n = 1)
sample_ids <- gsub('^"|"$', '', strsplit(hdr, " ")[[1]])
counts_dt <- fread(counts_file, sep = " ", header = FALSE, skip = 1,
                   data.table = FALSE, quote = '"')
gene_ids <- counts_dt[, 1]
counts_bulk <- as.matrix(counts_dt[, -1])
storage.mode(counts_bulk) <- "integer"
rownames(counts_bulk) <- gene_ids
colnames(counts_bulk) <- sample_ids

# Ensembl → Symbol
ens_clean <- sub("\\..*$", "", rownames(counts_bulk))
sym_map <- mapIds(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
new_rn <- ifelse(is.na(sym_map), ens_clean, sym_map)
v <- matrixStats::rowVars(counts_bulk)
ord <- order(v, decreasing = TRUE)
counts_bulk <- counts_bulk[ord, ]; new_rn <- new_rn[ord]
counts_bulk <- counts_bulk[!duplicated(new_rn), ]
rownames(counts_bulk) <- new_rn[!duplicated(new_rn)]
cat(sprintf("   Bulk counts: %d genes × %d samples\n",
            nrow(counts_bulk), ncol(counts_bulk)))

# Align with QC-passed samples (using same set as MuSiC analysis)
prop_music <- read.csv("results/04_MuSiC_cell_fractions.csv")
common_samples <- intersect(colnames(counts_bulk), prop_music$sample)
counts_bulk <- counts_bulk[, common_samples]
cat(sprintf("   QC-aligned: %d samples\n", ncol(counts_bulk)))

# Bulk ExpressionSet
bulk_eset <- ExpressionSet(assayData = counts_bulk)

# ---- 4. Bisque ReferenceBasedDecomposition ----
cat(">>> [3/6] Bisque ReferenceBasedDecomposition (수분 소요)\n")
# use.overlap = FALSE 필수 — Smillie donor가 GSE193677에 없음
# old.cpm = TRUE: scale to CPM before deconvolution

bisque_cache <- "rds/bisque_deconv.rds"
if (file.exists(bisque_cache)) {
  cat("   Loading cached Bisque result\n")
  bisque_res <- readRDS(bisque_cache)
} else {
  bisque_res <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset = bulk_eset,
    sc.eset = sc_eset,
    markers = NULL,        # Bisque가 자동으로 cell-type marker 선택
    cell.types = "cellType",
    subject.names = "SubjectName",
    use.overlap = FALSE,
    verbose = TRUE
  )
  saveRDS(bisque_res, bisque_cache)
}

# Bisque output: bisque_res$bulk.props — cell type × sample matrix
prop_bisque <- t(bisque_res$bulk.props)  # → sample × cellType
prop_bisque <- as.data.frame(prop_bisque)
prop_bisque$sample <- rownames(prop_bisque)

write.csv(prop_bisque, "results/17_Bisque_cell_fractions.csv", row.names = FALSE)
cat("Bisque cell fraction summary:\n")
print(summary(prop_bisque[, setdiff(colnames(prop_bisque), "sample")]))

# ---- 5. Cross-validation: MuSiC vs Bisque ----
cat(">>> [4/6] MuSiC vs Bisque cross-validation\n")

# 매칭
prop_music_e <- prop_music %>%
  dplyr::select(sample, Epithelial_music = Epithelial,
                Immune_music    = Immune)
prop_bisque_e <- prop_bisque %>%
  dplyr::select(sample, Epithelial_bisque = Epithelial,
                Immune_bisque    = Immune)
merged <- inner_join(prop_music_e, prop_bisque_e, by = "sample")
cat(sprintf("   Matched samples: %d\n", nrow(merged)))

# Spearman r (epithelial 우선, 사전 등록 threshold ≥ 0.7)
epi_spearman <- cor(merged$Epithelial_music, merged$Epithelial_bisque,
                    method = "spearman")
epi_pearson <- cor(merged$Epithelial_music, merged$Epithelial_bisque,
                   method = "pearson")
imm_spearman <- cor(merged$Immune_music, merged$Immune_bisque,
                    method = "spearman")
imm_pearson <- cor(merged$Immune_music, merged$Immune_bisque,
                   method = "pearson")

cat(sprintf("\n=== MuSiC ↔ Bisque agreement ===\n"))
cat(sprintf("  Epithelial fraction: Spearman r = %.4f, Pearson r = %.4f\n",
            epi_spearman, epi_pearson))
cat(sprintf("  Immune fraction:     Spearman r = %.4f, Pearson r = %.4f\n",
            imm_spearman, imm_pearson))
cat(sprintf("  Pre-registered threshold r ≥ 0.7 met? Epithelial: %s, Immune: %s\n",
            ifelse(epi_spearman >= 0.7, "YES ✓", "NO ✗"),
            ifelse(imm_spearman >= 0.7, "YES ✓", "NO ✗")))

# 결과 요약 저장
agreement_summary <- data.frame(
  cell_type = c("Epithelial", "Immune"),
  spearman_r = c(epi_spearman, imm_spearman),
  pearson_r = c(epi_pearson, imm_pearson),
  threshold_0.7_met = c(epi_spearman >= 0.7, imm_spearman >= 0.7),
  n_samples = nrow(merged)
)
write.csv(agreement_summary, "results/18_MuSiC_Bisque_agreement.csv",
          row.names = FALSE)

# ---- 6. Visualization — Supplementary Fig 5 ----
cat(">>> [5/6] Build Supplementary Fig 5\n")

# Epithelial scatter
p_epi <- ggplot(merged, aes(Epithelial_music, Epithelial_bisque)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(alpha = 0.3, size = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.6) +
  labs(x = "MuSiC Epithelial fraction",
       y = "Bisque Epithelial fraction",
       title = "MuSiC vs Bisque — Epithelial fraction agreement",
       subtitle = sprintf("Spearman r = %.4f, Pearson r = %.4f, n = %d samples",
                          epi_spearman, epi_pearson, nrow(merged))) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Immune scatter
p_imm <- ggplot(merged, aes(Immune_music, Immune_bisque)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(alpha = 0.3, size = 0.6, color = "firebrick") +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", linewidth = 0.6) +
  labs(x = "MuSiC Immune fraction",
       y = "Bisque Immune fraction",
       title = "MuSiC vs Bisque — Immune fraction agreement",
       subtitle = sprintf("Spearman r = %.4f, Pearson r = %.4f, n = %d samples",
                          imm_spearman, imm_pearson, nrow(merged))) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Combined panel
library(patchwork)
sfig5 <- p_epi + p_imm + plot_annotation(
  title = "Supplementary Figure 5 — Cross-validation of MuSiC cell fractions with Bisque",
  subtitle = "Same Smillie GSE116222 reference, GSE193677 bulk mixture (n=2,483); Pre-registered Spearman r ≥ 0.7 threshold"
)
ggsave("figures/21_supp5_MuSiC_vs_Bisque.png", sfig5,
       width = 12, height = 6, dpi = 150)

# Bland-Altman style plot (difference vs mean)
merged$mean_epi <- (merged$Epithelial_music + merged$Epithelial_bisque) / 2
merged$diff_epi <- merged$Epithelial_bisque - merged$Epithelial_music
mean_diff <- mean(merged$diff_epi)
sd_diff <- sd(merged$diff_epi)
p_ba <- ggplot(merged, aes(mean_epi, diff_epi)) +
  geom_point(alpha = 0.3, size = 0.6, color = "steelblue") +
  geom_hline(yintercept = mean_diff, color = "darkred") +
  geom_hline(yintercept = mean_diff + 1.96 * sd_diff,
             linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = mean_diff - 1.96 * sd_diff,
             linetype = "dashed", color = "darkred") +
  labs(x = "Mean Epithelial fraction (MuSiC + Bisque) / 2",
       y = "Difference (Bisque − MuSiC)",
       title = "Bland-Altman: MuSiC vs Bisque Epithelial fraction",
       subtitle = sprintf("Mean diff = %.3f, ±1.96 SD = %.3f", mean_diff, 1.96 * sd_diff)) +
  theme_bw(base_size = 10)
ggsave("figures/22_BlandAltman_MuSiC_Bisque.png", p_ba,
       width = 7, height = 5, dpi = 150)

cat("\n=== DONE ===\n")
cat("  results/17_Bisque_cell_fractions.csv\n")
cat("  results/18_MuSiC_Bisque_agreement.csv ★\n")
cat("  figures/21_supp5_MuSiC_vs_Bisque.png ★ Supplementary Fig 5\n")
cat("  figures/22_BlandAltman_MuSiC_Bisque.png\n")
cat(sprintf("\nKey result: Epithelial Spearman r = %.4f (threshold ≥ 0.7: %s)\n",
            epi_spearman, ifelse(epi_spearman >= 0.7, "PASS", "FAIL")))
