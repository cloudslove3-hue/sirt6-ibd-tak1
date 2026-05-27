# ==============================================================================
# MuSiC deconvolution + partial correlation
# Smillie scRNA reference로 GSE193677 bulk의 epithelial fraction 추정
# → SIRT6 ~ NLRP3 | epithelial_fraction 의 partial Spearman 계산
# → 1차 paper Fig 4 (결정적 figure)
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(MuSiC); library(Biobase); library(SingleCellExperiment)
  library(Matrix); library(dplyr); library(tidyr); library(ggplot2)
  library(ggpubr); library(pheatmap); library(tibble); library(stringr)
})
set.seed(42)

# ==============================================================================
# 1. Load scRNA reference (Smillie, broad cell types)
# ==============================================================================
cat(">>> [1/6] Load scRNA reference\n")
sce_ref <- readRDS("rds/smillie_sce_for_music.rds")
cat(sprintf("   reference: %d cells x %d genes\n",
            ncol(sce_ref), nrow(sce_ref)))
cat("   cellType distribution:\n"); print(table(colData(sce_ref)$cellType))
cat("   subject distribution:\n"); print(table(colData(sce_ref)$subjectID))

# ==============================================================================
# 2. Load bulk expression (GSE193677, vst 결과)
# ==============================================================================
cat(">>> [2/6] Load bulk GSE193677\n")
expr_bulk <- readRDS("../pilot_analysis/rds/expr_vst.rds")  # gene × sample, vst
pdata_bulk <- readRDS("../pilot_analysis/rds/pdata_final.rds")

# MuSiC requires raw counts (not vst). 다시 raw counts 로드 필요.
# Pilot 분석에서 raw count는 따로 저장 안 함 — 다시 로드
cat("   Note: MuSiC requires raw counts. Re-loading from GSE193677 counts file...\n")

library(data.table)
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

# Ensembl → Symbol (pilot과 동일 매핑)
library(org.Hs.eg.db); library(AnnotationDbi)
ens_clean <- sub("\\..*$", "", rownames(counts_bulk))
sym_map <- mapIds(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
new_rn <- ifelse(is.na(sym_map), ens_clean, sym_map)
# variance-based dedup
library(matrixStats)
v <- matrixStats::rowVars(counts_bulk)
ord <- order(v, decreasing = TRUE)
counts_bulk <- counts_bulk[ord, ]; new_rn <- new_rn[ord]
keep_dedup <- !duplicated(new_rn)
counts_bulk <- counts_bulk[keep_dedup, ]
rownames(counts_bulk) <- new_rn[keep_dedup]

# QC pass samples only (vst 매트릭스의 컬럼과 일치)
common_samples <- intersect(colnames(counts_bulk), colnames(expr_bulk))
counts_bulk <- counts_bulk[, common_samples]
cat(sprintf("   bulk count matrix (after QC alignment): %d genes x %d samples\n",
            nrow(counts_bulk), ncol(counts_bulk)))

# ==============================================================================
# 3. Common gene 교집합 + MuSiC ExpressionSet 구성
# ==============================================================================
cat(">>> [3/6] Build MuSiC inputs\n")
common_genes <- intersect(rownames(counts_bulk), rownames(sce_ref))
cat(sprintf("   common genes: %d\n", length(common_genes)))

# Bulk: ExpressionSet
pheno_bulk <- AnnotatedDataFrame(
  data.frame(SampleID = colnames(counts_bulk), row.names = colnames(counts_bulk))
)
eset_bulk <- ExpressionSet(
  assayData = counts_bulk[common_genes, ],
  phenoData = pheno_bulk
)

# Reference: filter to common genes
sce_ref_f <- sce_ref[common_genes, ]

# ==============================================================================
# 4. MuSiC deconvolution (cached)
# ==============================================================================
cache_music <- "rds/music_deconv.rds"
if (file.exists(cache_music)) {
  cat(">>> [4/6] Loading cached MuSiC deconvolution\n")
  deconv <- readRDS(cache_music)
  prop_est <- as.data.frame(deconv$Est.prop.weighted)
  prop_est$sample <- rownames(prop_est)
} else {
  cat(">>> [4/6] MuSiC deconvolution (수 분 소요)\n")
  deconv <- music_prop(
    bulk.mtx    = exprs(eset_bulk),
    sc.sce      = sce_ref_f,
    clusters    = "cellType",
    samples     = "subjectID",
    select.ct   = unique(colData(sce_ref_f)$cellType),
    verbose     = TRUE
  )
  prop_est <- as.data.frame(deconv$Est.prop.weighted)
  prop_est$sample <- rownames(prop_est)
  write.csv(prop_est, "results/04_MuSiC_cell_fractions.csv", row.names = FALSE)
  saveRDS(deconv, cache_music)
}
cat("   Cell fraction summary:\n")
print(summary(prop_est[, setdiff(colnames(prop_est), "sample")]))

# Visualize fractions by disease
prop_long <- prop_est %>%
  pivot_longer(-sample, names_to = "cellType", values_to = "fraction") %>%
  left_join(pdata_bulk %>% dplyr::select(biopsy_id, disease, inflammation, location),
            by = c("sample" = "biopsy_id"))

p_frac <- ggplot(prop_long, aes(disease, fraction, fill = inflammation)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~cellType, scales = "free_y") +
  labs(title = "MuSiC cell fractions by disease × inflammation") +
  theme_bw(base_size = 11)
ggsave("figures/11_cell_fractions_boxplot.png", p_frac,
       width = 11, height = 7, dpi = 130)

# ==============================================================================
# 5. Partial correlation: SIRT6 ~ NLRP3 | epithelial_fraction
# ==============================================================================
cat(">>> [5/6] Partial correlation analysis\n")

# Merge expression + fraction
expr_vst <- expr_bulk[, common_samples]
sirt6_vec <- expr_vst["SIRT6", ]
epi_frac <- prop_est$Epithelial[match(common_samples, prop_est$sample)]
if (all(is.na(epi_frac))) {
  # cellType label might differ (e.g., "epithelial" lowercase)
  candidate_epi <- grep("epi", colnames(prop_est), ignore.case = TRUE, value = TRUE)
  epi_frac <- prop_est[[candidate_epi[1]]][match(common_samples, prop_est$sample)]
  cat("   Using epithelial column:", candidate_epi[1], "\n")
}
imm_frac <- {
  col_imm <- grep("imm", colnames(prop_est), ignore.case = TRUE, value = TRUE)
  if (length(col_imm)) prop_est[[col_imm[1]]][match(common_samples, prop_est$sample)]
  else rep(NA, length(common_samples))
}

# Function: partial Spearman via rank-residual regression
partial_spearman <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  res_x <- resid(lm(rx ~ rz))
  res_y <- resid(lm(ry ~ rz))
  ct <- cor.test(res_x, res_y, method = "pearson")
  c(rho = unname(ct$estimate), p = ct$p.value)
}

targets <- c("NLRP3", "CASP1", "IL1B", "PYCARD", "MAP3K7", "RELA", "FOXC1")
targets <- intersect(targets, rownames(expr_vst))

# Subgroups (consistent with pilot)
subgroups <- list(
  All_samples    = common_samples,
  Control        = common_samples[pdata_bulk$disease[match(common_samples, pdata_bulk$biopsy_id)] == "Control"],
  UC_inflamed    = common_samples[pdata_bulk$disease[match(common_samples, pdata_bulk$biopsy_id)] == "UC" &
                                  pdata_bulk$inflammation[match(common_samples, pdata_bulk$biopsy_id)] == "I"],
  CD_inflamed    = common_samples[pdata_bulk$disease[match(common_samples, pdata_bulk$biopsy_id)] == "CD" &
                                  pdata_bulk$inflammation[match(common_samples, pdata_bulk$biopsy_id)] == "I"]
)

# Compute marginal vs partial for each subgroup × target
pcor_res <- list()
for (gn in names(subgroups)) {
  samps <- subgroups[[gn]]
  if (length(samps) < 20) next
  idx <- match(samps, common_samples)
  s_sub <- sirt6_vec[idx]
  epi_sub <- epi_frac[idx]
  for (tg in targets) {
    t_sub <- expr_vst[tg, idx]
    marg <- cor.test(s_sub, t_sub, method = "spearman", exact = FALSE)
    par <- partial_spearman(s_sub, t_sub, epi_sub)
    pcor_res[[length(pcor_res) + 1]] <- data.frame(
      subset = gn, target = tg, n = length(samps),
      rho_marginal = unname(marg$estimate), p_marginal = marg$p.value,
      rho_partial  = par["rho"], p_partial = par["p"]
    )
  }
}
pcor_df <- bind_rows(pcor_res) %>%
  group_by(subset) %>%
  mutate(q_marginal = p.adjust(p_marginal, "BH"),
         q_partial  = p.adjust(p_partial,  "BH")) %>%
  ungroup() %>%
  mutate(attenuation = abs(rho_marginal) - abs(rho_partial))

write.csv(pcor_df, "results/05_partial_correlation.csv", row.names = FALSE)
cat("   Marginal vs partial correlation:\n"); print(pcor_df)

# ==============================================================================
# 6. Fig 4 prototype — attenuation visualization
# ==============================================================================
cat(">>> [6/6] Fig 4 prototype\n")

# Heatmap: rho_marginal vs rho_partial side-by-side
pdf_long <- pcor_df %>%
  dplyr::select(subset, target, rho_marginal, rho_partial) %>%
  pivot_longer(c(rho_marginal, rho_partial), names_to = "type", values_to = "rho")

p_compare <- ggplot(pdf_long, aes(target, rho, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~subset, ncol = 2) +
  scale_fill_manual(values = c(rho_marginal = "#B2182B", rho_partial = "#2166AC"),
                    labels = c("Marginal", "Partial | epithelial_fraction")) +
  labs(y = "Spearman ρ", x = NULL, fill = NULL,
       title = "SIRT6 ~ target gene: marginal vs cell-composition-adjusted partial correlation",
       subtitle = "If partial rho ≈ 0 → marginal driven by cell composition. If preserved → independent signal.") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figures/12_marginal_vs_partial.png", p_compare,
       width = 12, height = 7, dpi = 130)

# Heatmap form
m1 <- pcor_df %>% dplyr::select(subset, target, rho_marginal) %>%
  pivot_wider(names_from = target, values_from = rho_marginal) %>%
  column_to_rownames("subset") %>% as.matrix()
m2 <- pcor_df %>% dplyr::select(subset, target, rho_partial) %>%
  pivot_wider(names_from = target, values_from = rho_partial) %>%
  column_to_rownames("subset") %>% as.matrix()

png("figures/13_heatmap_marginal.png", 800, 400, res = 130)
print(pheatmap(m1, color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
               breaks = seq(-0.6, 0.6, length.out=101),
               main = "Marginal ρ", display_numbers = TRUE, fontsize_number = 9,
               cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()
png("figures/14_heatmap_partial.png", 800, 400, res = 130)
print(pheatmap(m2, color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
               breaks = seq(-0.6, 0.6, length.out=101),
               main = "Partial ρ | epithelial_fraction",
               display_numbers = TRUE, fontsize_number = 9,
               cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

# Interpretation summary
cat("\n=== INTERPRETATION ===\n")
cat("For each (subset, target):\n")
cat("  - If |rho_partial| < 0.5 × |rho_marginal|  → SIRT6 ~ target는 주로 cell composition으로 설명\n")
cat("  - If |rho_partial| ≈ |rho_marginal|         → cell composition과 무관한 신호 (직접 조절 가능성)\n")
cat("  - 'attenuation' 컬럼이 partial로 감소한 크기\n\n")

summary_per_target <- pcor_df %>%
  filter(subset %in% c("UC_inflamed", "CD_inflamed")) %>%
  group_by(target) %>%
  summarise(
    avg_marg = mean(abs(rho_marginal)),
    avg_part = mean(abs(rho_partial)),
    pct_retained = avg_part / avg_marg * 100
  ) %>%
  arrange(desc(avg_marg))
print(summary_per_target)

cat("\n=== DONE ===\n")
cat("  - results/04_MuSiC_cell_fractions.csv\n")
cat("  - results/05_partial_correlation.csv ★\n")
cat("  - figures/11~14.png (cell fractions, marginal vs partial)\n")
cat("  - figures/12_marginal_vs_partial.png ★ Fig 4 prototype\n")
