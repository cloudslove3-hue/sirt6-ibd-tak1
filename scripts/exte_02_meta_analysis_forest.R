# ==============================================================================
# Fig 5C/5D 보강 — Forest plot + random-effects meta-analysis
# 사전 등록 (S1 §4-4): "≥60% direction-consistency in 28 discovery-validation
# pairs" 충족 후 본 meta-analysis 진행
# Method: Fisher z transformation → escalc → rma (random-effects) → forest
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/external_validation")

suppressPackageStartupMessages({
  library(metafor); library(dplyr); library(tidyr); library(ggplot2)
  library(tibble); library(stringr); library(ggrepel); library(scales)
})

# ==============================================================================
# 1. Load discovery + validation correlation tables
# ==============================================================================
cat(">>> [1/5] Load correlation tables\n")
disc_full <- read.csv("../pilot_analysis/results/01_sirt6_correlation_table.csv")
val_full  <- read.csv("results/01_sirt6_correlation.csv")

# 4 subgroup만 비교 가능 (val에서 subgroup은 All/Control/UC/CD)
# Discovery: All_samples, Control, UC_all (=UC), CD_all (=CD) 선택
disc <- disc_full %>%
  filter(subset %in% c("All_samples", "Control", "UC_all", "CD_all")) %>%
  mutate(subset_clean = case_when(
    subset == "UC_all" ~ "UC",
    subset == "CD_all" ~ "CD",
    TRUE ~ subset)) %>%
  dplyr::select(subset = subset_clean, n, target_label, target_gene, rho)
disc$study <- "Discovery (GSE193677)"

val <- val_full %>%
  filter(subset %in% c("All_samples", "Control", "UC", "CD")) %>%
  dplyr::select(subset, n, target_label, target_gene, rho)
val$study <- "Validation (GSE235236)"

all_es <- bind_rows(disc, val)
cat(sprintf("   total effect sizes: %d (subgroup × target × study)\n",
            nrow(all_es)))

# ==============================================================================
# 2. Fisher z transformation + variance
# ==============================================================================
cat(">>> [2/5] Fisher z + variance\n")
# escalc(measure='ZCOR', ri=rho, ni=n) — Spearman ρ를 Pearson-like Fisher z로 근사
# Spearman의 Fisher z 변환은 보수적 근사이며, 표준 sample variance: 1/(n-3)
es <- escalc(measure = "ZCOR", ri = all_es$rho, ni = all_es$n,
             data = all_es)
es$study <- all_es$study
es$subset <- all_es$subset
es$target_label <- all_es$target_label
# back-transform yi (Fisher z) → ρ
es$rho_orig <- (exp(2 * es$yi) - 1) / (exp(2 * es$yi) + 1)
es$rho_lower <- (exp(2 * (es$yi - 1.96 * sqrt(es$vi))) - 1) /
                (exp(2 * (es$yi - 1.96 * sqrt(es$vi))) + 1)
es$rho_upper <- (exp(2 * (es$yi + 1.96 * sqrt(es$vi))) - 1) /
                (exp(2 * (es$yi + 1.96 * sqrt(es$vi))) + 1)

write.csv(es %>% dplyr::select(study, subset, target_label, n,
                                rho_orig, rho_lower, rho_upper, yi, vi),
          "results/05_meta_effect_sizes.csv", row.names = FALSE)

# ==============================================================================
# 3. Per-target random-effects meta-analysis
# ==============================================================================
cat(">>> [3/5] Per-target random-effects meta-analysis (REML)\n")
targets <- unique(es$target_label)
meta_results <- list()
for (tg in targets) {
  sub <- es %>% filter(target_label == tg)
  if (nrow(sub) < 3) next
  m <- rma(yi = sub$yi, vi = sub$vi, method = "REML")
  pooled_z <- as.numeric(m$beta)
  pooled_rho <- (exp(2 * pooled_z) - 1) / (exp(2 * pooled_z) + 1)
  pooled_ci_lo <- (exp(2 * m$ci.lb) - 1) / (exp(2 * m$ci.lb) + 1)
  pooled_ci_hi <- (exp(2 * m$ci.ub) - 1) / (exp(2 * m$ci.ub) + 1)
  meta_results[[tg]] <- data.frame(
    target = tg,
    k = nrow(sub),
    pooled_z = pooled_z,
    pooled_rho = pooled_rho,
    pooled_lower = pooled_ci_lo,
    pooled_upper = pooled_ci_hi,
    p_pooled = m$pval,
    I2 = m$I2,
    tau2 = m$tau2,
    Q_p = m$QEp
  )
}
meta_df <- bind_rows(meta_results) %>%
  arrange(pooled_rho)
write.csv(meta_df, "results/06_meta_analysis_pooled.csv", row.names = FALSE)
cat("Pooled per-target (random-effects):\n"); print(meta_df)

# ==============================================================================
# 4. Forest plot — 8 row per target × 7 target panels
# ==============================================================================
cat(">>> [4/5] Build forest plot\n")

# Order subset for consistent display
es$subset <- factor(es$subset, levels = c("All_samples", "Control", "UC", "CD"))
es$study <- factor(es$study, levels = c("Discovery (GSE193677)", "Validation (GSE235236)"))
es$row_id <- paste(es$subset, es$study, sep = " | ")
target_order <- meta_df$target  # sorted by pooled_rho ascending
es$target_label <- factor(es$target_label, levels = target_order)

p_forest <- ggplot(es, aes(y = row_id, x = rho_orig,
                            color = study)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = rho_lower, xmax = rho_upper), height = 0.3) +
  geom_point(aes(size = n)) +
  geom_vline(data = meta_df %>%
               rename(target_label = target),
             aes(xintercept = pooled_rho),
             linetype = "dotted", color = "firebrick", linewidth = 0.6) +
  geom_rect(data = meta_df %>% rename(target_label = target),
            aes(xmin = pooled_lower, xmax = pooled_upper,
                ymin = -Inf, ymax = Inf),
            fill = "firebrick", alpha = 0.08, inherit.aes = FALSE) +
  facet_wrap(~target_label, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Discovery (GSE193677)" = "#1f78b4",
                                 "Validation (GSE235236)" = "#e31a1c")) +
  scale_size_continuous(range = c(1.5, 5), trans = "log10",
                        breaks = c(10, 100, 1000)) +
  labs(x = "Spearman ρ (95% CI, Fisher z)",
       y = NULL,
       title = "SIRT6 vs target gene — random-effects meta-analysis",
       subtitle = "Each row: subgroup × study. Red dotted line: pooled ρ (REML).",
       color = NULL, size = "n samples") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("figures/03_meta_forest.png", p_forest,
       width = 14, height = 9, dpi = 150)

# ==============================================================================
# 5. Pooled summary plot (Fig 5D) — pooled ρ per target with I²
# ==============================================================================
cat(">>> [5/5] Build Fig 5D — pooled summary\n")
meta_plot <- meta_df %>%
  mutate(target = factor(target, levels = rev(target))) %>%
  mutate(label_I2 = sprintf("I² = %.0f%%", I2),
         label_n = sprintf("k = %d", k))

p_pool <- ggplot(meta_plot, aes(y = target, x = pooled_rho)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = pooled_lower, xmax = pooled_upper),
                 height = 0.3, color = "black", linewidth = 0.6) +
  geom_point(size = 4, color = "firebrick") +
  geom_text(aes(label = label_I2, x = pooled_upper + 0.05),
            hjust = 0, size = 3.2, color = "gray30") +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.25)) +
  labs(x = "Pooled Spearman ρ (random-effects, 95% CI)",
       y = NULL,
       title = "Cross-cohort pooled effect sizes",
       subtitle = "Random-effects REML; k = subgroup × study comparisons per target") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("figures/04_meta_pooled_summary.png", p_pool,
       width = 9, height = 5, dpi = 150)

# ==============================================================================
# 6. Subgroup-specific pooled (Discovery vs Validation only — sensitivity)
# ==============================================================================
cat(">>> [bonus] Per-study pooled for TAK1 / NLRP3\n")
for (tg in c("TAK1", "NLRP3")) {
  sub <- es %>% filter(target_label == tg)
  m_disc <- rma(yi = sub$yi[sub$study == "Discovery (GSE193677)"],
                vi = sub$vi[sub$study == "Discovery (GSE193677)"],
                method = "REML")
  m_val <- rma(yi = sub$yi[sub$study == "Validation (GSE235236)"],
               vi = sub$vi[sub$study == "Validation (GSE235236)"],
               method = "REML")
  pooled_d <- as.numeric(m_disc$beta); pooled_v <- as.numeric(m_val$beta)
  rho_d <- (exp(2*pooled_d)-1)/(exp(2*pooled_d)+1)
  rho_v <- (exp(2*pooled_v)-1)/(exp(2*pooled_v)+1)
  cat(sprintf("  %s: Discovery pooled ρ=%.3f (I²=%.0f%%), Validation pooled ρ=%.3f (I²=%.0f%%)\n",
              tg, rho_d, m_disc$I2, rho_v, m_val$I2))
}

cat("\n=== DONE ===\n")
cat("Outputs:\n")
cat("  results/05_meta_effect_sizes.csv — 56 effect sizes (subgroup × target × study)\n")
cat("  results/06_meta_analysis_pooled.csv — per-target pooled ρ + I²\n")
cat("  figures/03_meta_forest.png ★ Fig 5C — 56 effect-size forest plot\n")
cat("  figures/04_meta_pooled_summary.png ★ Fig 5D — pooled per-target with I²\n")
