# ==============================================================================
# Supplementary — ppcor sensitivity check
# 본 분석의 rank-residual partial Spearman vs ppcor 패키지 표준 partial Spearman
# 일치도 확인 (Methods §4d, manuscript v2)
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(ppcor); library(dplyr); library(tidyr); library(ggplot2)
  library(tibble)
})

# Load expr_vst + pdata + MuSiC fractions
expr <- readRDS("../pilot_analysis/rds/expr_vst.rds")
pdata <- readRDS("../pilot_analysis/rds/pdata_final.rds")
prop_est <- read.csv("results/04_MuSiC_cell_fractions.csv")

common <- intersect(colnames(expr), prop_est$sample)
expr <- expr[, common]
pdata_sub <- pdata[match(common, pdata$biopsy_id), ]
epi_col <- grep("epi", colnames(prop_est), ignore.case = TRUE, value = TRUE)[1]
epi_frac <- prop_est[[epi_col]][match(common, prop_est$sample)]

targets <- c("NLRP3", "CASP1", "IL1B", "PYCARD", "MAP3K7", "RELA", "FOXC1")
targets <- intersect(targets, rownames(expr))

# Subgroups
subgroups <- list(
  All_samples = common,
  Control     = common[pdata_sub$disease == "Control"],
  UC_inflamed = common[pdata_sub$disease == "UC" & pdata_sub$inflammation == "I"],
  CD_inflamed = common[pdata_sub$disease == "CD" & pdata_sub$inflammation == "I"]
)

# Rank-residual partial (본 분석 method)
partial_rankres <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  cor(resid(lm(rx ~ rz)), resid(lm(ry ~ rz)), method = "pearson")
}

cat(">>> Compute both methods for each (subgroup, target)\n")
res <- list()
for (gn in names(subgroups)) {
  samps <- subgroups[[gn]]
  if (length(samps) < 20) next
  idx <- match(samps, common)
  s <- expr["SIRT6", idx]; e <- epi_frac[idx]
  for (tg in targets) {
    t <- expr[tg, idx]
    rho_rankres <- partial_rankres(s, t, e)
    # ppcor: pcor.test with Spearman method
    pp <- ppcor::pcor.test(s, t, e, method = "spearman")
    res[[length(res)+1]] <- data.frame(
      subset = gn, target = tg, n = length(samps),
      rho_rankres = rho_rankres,
      rho_ppcor = pp$estimate,
      diff = rho_rankres - pp$estimate
    )
  }
}
res_df <- bind_rows(res)
write.csv(res_df, "results/14_ppcor_sensitivity.csv", row.names = FALSE)

# Agreement statistics
agreement_pearson <- cor(res_df$rho_rankres, res_df$rho_ppcor, method = "pearson")
agreement_spearman <- cor(res_df$rho_rankres, res_df$rho_ppcor, method = "spearman")
max_abs_diff <- max(abs(res_df$diff))
cat(sprintf("\nAgreement Pearson r = %.4f\n", agreement_pearson))
cat(sprintf("Agreement Spearman r = %.4f\n", agreement_spearman))
cat(sprintf("Max absolute difference = %.4f\n", max_abs_diff))

# Scatter plot
p <- ggplot(res_df, aes(x = rho_ppcor, y = rho_rankres, color = subset)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = target), size = 3, max.overlaps = 30) +
  labs(x = "ppcor partial Spearman ρ",
       y = "rank-residual partial Spearman ρ",
       title = "Partial correlation method sensitivity check",
       subtitle = sprintf("Pearson r = %.4f, Spearman r = %.4f, max |diff| = %.4f",
                          agreement_pearson, agreement_spearman, max_abs_diff)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/19_ppcor_sensitivity.png", p, width = 8, height = 6, dpi = 130)

cat("\n=== DONE ===\n")
cat("  results/14_ppcor_sensitivity.csv\n")
cat("  figures/19_ppcor_sensitivity.png ★ Supplementary Fig SX\n")
