# ==============================================================================
# Bootstrap 95% CI for partial Spearman (per subgroup × target)
# 사용자 권고 (b) — ppcor sensitivity의 algebraic equivalence를 대체할
# 실질적 robustness check
# Method: 각 subgroup 안에서 sample bootstrap (n_boot=1000) → 2.5%/97.5% quantile
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
})

set.seed(42)
N_BOOT <- 1000

# Load
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
subgroups <- list(
  All_samples = common,
  Control     = common[pdata_sub$disease == "Control"],
  UC_inflamed = common[pdata_sub$disease == "UC" & pdata_sub$inflammation == "I"],
  CD_inflamed = common[pdata_sub$disease == "CD" & pdata_sub$inflammation == "I"]
)

partial_rankres <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  cor(resid(lm(rx ~ rz)), resid(lm(ry ~ rz)), method = "pearson")
}

cat(sprintf(">>> Bootstrap: N=%d resamples per (subgroup × target)\n", N_BOOT))
cat(sprintf("    targets: %d, subgroups: %d → %d combinations\n",
            length(targets), length(subgroups), length(targets)*length(subgroups)))

results <- list()
for (gn in names(subgroups)) {
  samps <- subgroups[[gn]]
  if (length(samps) < 20) next
  idx <- match(samps, common)
  s_obs <- expr["SIRT6", idx]
  e_obs <- epi_frac[idx]
  for (tg in targets) {
    t_obs <- expr[tg, idx]
    rho_point <- partial_rankres(s_obs, t_obs, e_obs)
    # Bootstrap
    boot_rhos <- replicate(N_BOOT, {
      bidx <- sample(seq_along(s_obs), replace = TRUE)
      partial_rankres(s_obs[bidx], t_obs[bidx], e_obs[bidx])
    })
    ci <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)
    results[[length(results)+1]] <- data.frame(
      subset = gn, target = tg, n = length(samps),
      rho_point = rho_point,
      ci_lower = unname(ci[1]),
      ci_upper = unname(ci[2]),
      boot_sd = sd(boot_rhos, na.rm = TRUE),
      boot_iqr = IQR(boot_rhos, na.rm = TRUE)
    )
  }
}
res_df <- bind_rows(results)
write.csv(res_df, "results/16_bootstrap_partial_CI.csv", row.names = FALSE)
cat("\nBootstrap CI table:\n")
print(as.data.frame(res_df))

# Visualization — forest plot per subgroup
res_df$target <- factor(res_df$target, levels = unique(res_df$target))
res_df$subset <- factor(res_df$subset,
                        levels = c("All_samples", "Control", "UC_inflamed", "CD_inflamed"))

p <- ggplot(res_df, aes(y = target, x = rho_point, color = subset)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_point(size = 3) +
  facet_wrap(~subset, ncol = 2) +
  scale_color_manual(values = c("All_samples"="black", "Control"="steelblue",
                                "UC_inflamed"="firebrick", "CD_inflamed"="darkorange")) +
  labs(x = "Partial Spearman ρ | epithelial_fraction (95% bootstrap CI, 1000 resamples)",
       y = NULL,
       title = "Bootstrap 95% CI for partial Spearman",
       subtitle = "Robustness check substituting algebraically-equivalent ppcor sensitivity") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave("figures/20_bootstrap_partial_CI.png", p, width = 11, height = 7, dpi = 150)

# 핵심 target의 CI 폭 보고
cat("\n=== Key target robustness ===\n")
key_summary <- res_df %>%
  filter(target %in% c("MAP3K7", "NLRP3", "RELA")) %>%
  mutate(ci_width = ci_upper - ci_lower) %>%
  dplyr::select(subset, target, rho_point, ci_lower, ci_upper, ci_width)
print(as.data.frame(key_summary))

cat("\n=== DONE ===\n")
cat("  results/16_bootstrap_partial_CI.csv\n")
cat("  figures/20_bootstrap_partial_CI.png ★ Supplementary Fig (replaces ppcor sensitivity)\n")
