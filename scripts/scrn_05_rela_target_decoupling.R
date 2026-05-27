# ==============================================================================
# RELA target gene 상관 — SIRT6 mRNA-activity decoupling 가설 검증
# 1차 paper Supplementary
#
# 가설: SIRT6와 RELA는 같은 epithelial program에 묶여 +ρ를 보이지만,
#       SIRT6는 RELA의 enzymatic transactivation을 억제 → RELA target gene
#       (IL6/TNF/CCL2/CXCL8)은 SIRT6와 *음의* 상관을 보일 것.
# 검증: subgroup별 SIRT6 vs RELA target genes 상관 + partial(epi 보정) 비교.
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(pheatmap)
  library(tibble)
})

# Inputs
expr_vst <- readRDS("../pilot_analysis/rds/expr_vst.rds")
pdata    <- readRDS("../pilot_analysis/rds/pdata_final.rds")
prop_est <- read.csv("results/04_MuSiC_cell_fractions.csv")

# RELA canonical target genes (NF-κB pro-inflammatory program)
rela_targets <- c("IL6", "TNF", "CCL2", "CXCL8", "IL8", "CXCL1", "CXCL2",
                  "ICAM1", "VCAM1", "NFKBIA", "TNFAIP3",   # 자가 조절 (negative feedback)
                  "CCL20", "PTGS2",  # COX2
                  "BCL2", "BIRC3", "XIAP")  # anti-apoptosis
rela_targets <- intersect(rela_targets, rownames(expr_vst))
cat("Available RELA targets:", paste(rela_targets, collapse=", "), "\n")

# Subgroups
common <- intersect(colnames(expr_vst), prop_est$sample)
pdata_sub <- pdata[match(common, pdata$biopsy_id), ]

# epi fraction lookup
epi_col <- grep("epi", colnames(prop_est), ignore.case = TRUE, value = TRUE)[1]
epi_frac <- prop_est[[epi_col]][match(common, prop_est$sample)]

subgroups <- list(
  All_samples    = common,
  Control        = common[pdata_sub$disease == "Control"],
  UC_inflamed    = common[pdata_sub$disease == "UC" & pdata_sub$inflammation == "I"],
  CD_inflamed    = common[pdata_sub$disease == "CD" & pdata_sub$inflammation == "I"]
)

partial_spearman <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  res_x <- resid(lm(rx ~ rz))
  res_y <- resid(lm(ry ~ rz))
  ct <- cor.test(res_x, res_y, method = "pearson")
  c(rho = unname(ct$estimate), p = ct$p.value)
}

# Compute SIRT6 vs each RELA target — marginal + partial
results <- list()
sirt6_all <- expr_vst["SIRT6", common]
for (gn in names(subgroups)) {
  samps <- subgroups[[gn]]
  if (length(samps) < 20) next
  idx <- match(samps, common)
  s <- sirt6_all[idx]; e <- epi_frac[idx]
  for (tg in rela_targets) {
    t <- expr_vst[tg, samps]
    marg <- cor.test(s, t, method = "spearman", exact = FALSE)
    par <- partial_spearman(s, t, e)
    results[[length(results) + 1]] <- data.frame(
      subset = gn, target = tg, n = length(samps),
      rho_marg = unname(marg$estimate), p_marg = marg$p.value,
      rho_part = par["rho"], p_part = par["p"]
    )
  }
}
res_df <- bind_rows(results) %>%
  group_by(subset) %>%
  mutate(q_marg = p.adjust(p_marg, "BH"),
         q_part = p.adjust(p_part, "BH")) %>% ungroup()

write.csv(res_df, "results/10_RELA_targets_partial.csv", row.names = FALSE)

# Heatmap
cat("\nMarginal SIRT6 vs RELA targets:\n")
m_marg <- res_df %>% select(subset, target, rho_marg) %>%
  pivot_wider(names_from = target, values_from = rho_marg) %>%
  column_to_rownames("subset") %>% as.matrix()
m_part <- res_df %>% select(subset, target, rho_part) %>%
  pivot_wider(names_from = target, values_from = rho_part) %>%
  column_to_rownames("subset") %>% as.matrix()

png("figures/15_RELA_targets_marginal.png", 1100, 400, res = 130)
print(pheatmap(m_marg,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks = seq(-0.5, 0.5, length.out = 101),
  main = "SIRT6 vs RELA target — marginal ρ",
  display_numbers = TRUE, fontsize_number = 8,
  cluster_rows = FALSE, cluster_cols = TRUE))
dev.off()

png("figures/16_RELA_targets_partial.png", 1100, 400, res = 130)
print(pheatmap(m_part,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks = seq(-0.5, 0.5, length.out = 101),
  main = "SIRT6 vs RELA target — partial ρ | epithelial_fraction",
  display_numbers = TRUE, fontsize_number = 8,
  cluster_rows = FALSE, cluster_cols = TRUE))
dev.off()

# Interpretation
cat("\nInterpretation key:\n")
cat("  - RELA(p65) 자체와는 ρ ≈ +0.25 (pilot 결과)\n")
cat("  - 그러나 RELA target gene이 SIRT6와 *음*의 상관 (또는 partial에서 음)이면\n")
cat("    → 'mRNA decoupling' 가설 지지: SIRT6 mRNA는 RELA mRNA와 동조하지만,\n")
cat("       enzymatic activity로 RELA transactivation을 억제\n")
cat("  - 만약 RELA target도 SIRT6와 양의 상관이면 → epithelial co-regulation 가설\n\n")

# 요약 통계
summary_df <- res_df %>%
  filter(subset %in% c("UC_inflamed", "CD_inflamed")) %>%
  group_by(target) %>%
  summarise(
    mean_marg = mean(rho_marg),
    mean_part = mean(rho_part),
    direction = case_when(
      mean_part < -0.1 ~ "Decoupled (SIRT6 suppresses target)",
      mean_part > 0.1  ~ "Co-regulated (epithelial program)",
      TRUE             ~ "Independent"
    )
  )
print(summary_df)
write.csv(summary_df, "results/11_RELA_targets_summary.csv", row.names = FALSE)
