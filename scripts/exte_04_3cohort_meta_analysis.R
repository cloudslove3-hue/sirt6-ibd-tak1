# ==============================================================================
# 3-cohort meta-analysis
# Cohorts: GSE193677 (discovery, n=2,483, RNA-Seq)
#          GSE235236 (val 1, n=56, RNA-Seq)
#          GSE75214  (val 2, n=194, microarray) — added by S1 Amendment 02
# Method: Random-effects REML, Fisher z, escalc(measure="ZCOR")
# 추가 sensitivity: RNA-Seq only subset (Discovery + GSE235236)
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/external_validation")

suppressPackageStartupMessages({
  library(metafor); library(dplyr); library(tidyr); library(ggplot2)
  library(tibble); library(stringr); library(ggrepel); library(scales); library(pheatmap)
})

# ==============================================================================
# 1. Load all three cohort correlations
# ==============================================================================
cat(">>> [1/6] Load 3 cohort correlation tables\n")
disc_full <- read.csv("../pilot_analysis/results/01_sirt6_correlation_table.csv")
val1_full <- read.csv("results/01_sirt6_correlation.csv")
val2_full <- read.csv("results_GSE75214/01_sirt6_correlation.csv")

disc <- disc_full %>%
  filter(subset %in% c("All_samples", "Control", "UC_all", "CD_all")) %>%
  mutate(subset = case_when(subset=="UC_all"~"UC", subset=="CD_all"~"CD", TRUE~subset)) %>%
  dplyr::select(subset, n, target_label, target_gene, rho)
disc$study <- "Discovery (GSE193677)"
disc$platform <- "RNA-Seq"

val1 <- val1_full %>%
  filter(subset %in% c("All_samples", "Control", "UC", "CD")) %>%
  dplyr::select(subset, n, target_label, target_gene, rho)
val1$study <- "Validation #1 (GSE235236)"
val1$platform <- "RNA-Seq"

val2 <- val2_full %>%
  filter(subset %in% c("All_samples", "Control", "UC", "CD")) %>%
  dplyr::select(subset, n, target_label, target_gene, rho)
val2$study <- "Validation #2 (GSE75214)"
val2$platform <- "Microarray"

all_es <- bind_rows(disc, val1, val2)
cat(sprintf("   total effect sizes: %d\n", nrow(all_es)))
cat("   per study:\n"); print(table(all_es$study))

# ==============================================================================
# 2. Fisher z + variance
# ==============================================================================
cat(">>> [2/6] Fisher z transformation\n")
es <- escalc(measure = "ZCOR", ri = all_es$rho, ni = all_es$n, data = all_es)
es$study <- all_es$study
es$subset <- all_es$subset
es$target_label <- all_es$target_label
es$platform <- all_es$platform
es$rho_orig <- (exp(2*es$yi)-1)/(exp(2*es$yi)+1)
es$rho_lower <- (exp(2*(es$yi-1.96*sqrt(es$vi)))-1)/(exp(2*(es$yi-1.96*sqrt(es$vi)))+1)
es$rho_upper <- (exp(2*(es$yi+1.96*sqrt(es$vi)))-1)/(exp(2*(es$yi+1.96*sqrt(es$vi)))+1)
write.csv(es %>% dplyr::select(study, platform, subset, target_label, n,
                                rho_orig, rho_lower, rho_upper, yi, vi),
          "results/08_3cohort_effect_sizes.csv", row.names = FALSE)

# ==============================================================================
# 3. Per-target meta-analysis — All studies pooled
# ==============================================================================
cat(">>> [3/6] All-3-cohort REML meta-analysis\n")
targets <- unique(es$target_label)
meta_all <- list()
for (tg in targets) {
  sub <- es %>% filter(target_label == tg)
  if (nrow(sub) < 3) next
  m <- rma(yi = sub$yi, vi = sub$vi, method = "REML")
  z <- as.numeric(m$beta)
  meta_all[[tg]] <- data.frame(
    target = tg, k = nrow(sub), pooled_z = z,
    pooled_rho   = (exp(2*z)-1)/(exp(2*z)+1),
    pooled_lower = (exp(2*m$ci.lb)-1)/(exp(2*m$ci.lb)+1),
    pooled_upper = (exp(2*m$ci.ub)-1)/(exp(2*m$ci.ub)+1),
    p_pooled = m$pval, I2 = m$I2, tau2 = m$tau2, Q_p = m$QEp,
    scope = "All 3 cohorts"
  )
}
meta_3 <- bind_rows(meta_all)

# ==============================================================================
# 4. Sensitivity: RNA-Seq only (Discovery + GSE235236)
# ==============================================================================
cat(">>> [4/6] Sensitivity: RNA-Seq-only pool\n")
es_rnaseq <- es %>% filter(platform == "RNA-Seq")
meta_rnaseq <- list()
for (tg in targets) {
  sub <- es_rnaseq %>% filter(target_label == tg)
  if (nrow(sub) < 3) next
  m <- rma(yi = sub$yi, vi = sub$vi, method = "REML")
  z <- as.numeric(m$beta)
  meta_rnaseq[[tg]] <- data.frame(
    target = tg, k = nrow(sub), pooled_z = z,
    pooled_rho   = (exp(2*z)-1)/(exp(2*z)+1),
    pooled_lower = (exp(2*m$ci.lb)-1)/(exp(2*m$ci.lb)+1),
    pooled_upper = (exp(2*m$ci.ub)-1)/(exp(2*m$ci.ub)+1),
    p_pooled = m$pval, I2 = m$I2, tau2 = m$tau2, Q_p = m$QEp,
    scope = "RNA-Seq only"
  )
}
meta_rna <- bind_rows(meta_rnaseq)

# ==============================================================================
# 5. Per-cohort: 4 subgroup × 1 study pooled (3 cohort separate)
# ==============================================================================
cat(">>> [5/6] Per-cohort pooled (4 subgroup → 1 estimate per cohort × target)\n")
meta_per_cohort <- list()
for (tg in targets) {
  for (st in unique(es$study)) {
    sub <- es %>% filter(target_label == tg, study == st)
    if (nrow(sub) < 2) next
    m <- rma(yi = sub$yi, vi = sub$vi, method = "REML")
    z <- as.numeric(m$beta)
    meta_per_cohort[[length(meta_per_cohort)+1]] <- data.frame(
      target = tg, study = st, k = nrow(sub), pooled_z = z,
      pooled_rho   = (exp(2*z)-1)/(exp(2*z)+1),
      pooled_lower = (exp(2*m$ci.lb)-1)/(exp(2*m$ci.lb)+1),
      pooled_upper = (exp(2*m$ci.ub)-1)/(exp(2*m$ci.ub)+1),
      I2 = m$I2
    )
  }
}
meta_per <- bind_rows(meta_per_cohort)
write.csv(meta_per, "results/09_per_cohort_pooled.csv", row.names = FALSE)

# Combine for final output
meta_combined <- bind_rows(meta_3, meta_rna) %>% arrange(target, scope)
write.csv(meta_combined, "results/10_3cohort_meta_pooled.csv", row.names = FALSE)
cat("\n=== Per-target pooled — All 3 vs RNA-Seq only ===\n")
print(as.data.frame(meta_combined))

# ==============================================================================
# 6. Visualizations — updated Fig 5C/5D
# ==============================================================================
cat(">>> [6/6] Updated figures\n")

# Order targets by All-3-cohort pooled rho
target_order <- meta_3 %>% arrange(pooled_rho) %>% pull(target)
es$target_label <- factor(es$target_label, levels = target_order)
es$subset <- factor(es$subset, levels = c("All_samples", "Control", "UC", "CD"))
es$study <- factor(es$study, levels = c("Discovery (GSE193677)",
                                          "Validation #1 (GSE235236)",
                                          "Validation #2 (GSE75214)"))

# Fig 5C: 84 effect-size forest plot (3 cohort × 4 subgroup × 7 target)
es$row_id <- paste(es$subset, sub(".* \\(", "(", es$study), sep = " | ")

p_forest <- ggplot(es, aes(y = row_id, x = rho_orig, color = study, shape = platform)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = rho_lower, xmax = rho_upper), height = 0.25) +
  geom_point(aes(size = n)) +
  geom_vline(data = meta_3 %>% rename(target_label = target),
             aes(xintercept = pooled_rho),
             linetype = "dotted", color = "firebrick", linewidth = 0.6) +
  facet_wrap(~target_label, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("Discovery (GSE193677)" = "#1f78b4",
                                 "Validation #1 (GSE235236)" = "#e31a1c",
                                 "Validation #2 (GSE75214)" = "#33a02c")) +
  scale_shape_manual(values = c("RNA-Seq" = 16, "Microarray" = 17)) +
  scale_size_continuous(range = c(1.2, 5), trans = "log10",
                        breaks = c(10, 100, 1000)) +
  labs(x = "Spearman ρ (95% CI)", y = NULL, color = NULL, shape = "Platform", size = "n",
       title = "3-cohort meta-analysis — SIRT6 vs target gene",
       subtitle = "Discovery (RNA-Seq, n=2,483) + GSE235236 (RNA-Seq, n=56) + GSE75214 (microarray, n=194). Red dotted = REML pooled.") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(face = "italic"),
        legend.position = "bottom", legend.box = "vertical",
        plot.title = element_text(face = "bold"))
ggsave("figures/06_3cohort_meta_forest.png", p_forest,
       width = 16, height = 10, dpi = 150)

# Fig 5D: Per-target pooled — All vs RNA-Seq only (Sensitivity)
meta_combined$target <- factor(meta_combined$target, levels = rev(target_order))
meta_combined$scope <- factor(meta_combined$scope, levels = c("RNA-Seq only", "All 3 cohorts"))
meta_combined$y_offset <- ifelse(meta_combined$scope == "All 3 cohorts", 0.18, -0.18)
meta_combined$y_pos <- as.numeric(meta_combined$target) + meta_combined$y_offset

p_pool <- ggplot(meta_combined, aes(y = y_pos, x = pooled_rho, color = scope)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = pooled_lower, xmax = pooled_upper), height = 0.1, linewidth = 0.7) +
  geom_point(size = 4, aes(shape = scope)) +
  geom_text(aes(label = sprintf("I²=%.0f%%, k=%d", I2, k),
                x = pmax(pooled_upper + 0.05, 0.5)),
            hjust = 0, size = 2.8, color = "gray30") +
  scale_y_continuous(breaks = seq_along(levels(meta_combined$target)),
                     labels = levels(meta_combined$target),
                     expand = expansion(add = 0.4)) +
  scale_x_continuous(limits = c(-1, 1.1), breaks = seq(-1, 1, 0.25)) +
  scale_color_manual(values = c("All 3 cohorts" = "black", "RNA-Seq only" = "#1f78b4")) +
  scale_shape_manual(values = c("All 3 cohorts" = 18, "RNA-Seq only" = 16)) +
  labs(x = "Pooled Spearman ρ (random-effects REML, 95% CI)", y = NULL,
       color = NULL, shape = NULL,
       title = "Pooled effect sizes — 3-cohort vs RNA-Seq-only sensitivity",
       subtitle = "GSE75214 (microarray) addition heterogeneity impact. I² shown for each scope.") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave("figures/07_3cohort_pooled_sensitivity.png", p_pool,
       width = 11, height = 6, dpi = 150)

# ==============================================================================
# 7. Direction consistency summary (Discovery vs 2 validations)
# ==============================================================================
cat("\n=== Direction-consistency summary ===\n")
# Discovery as reference
disc_ref <- disc %>%
  dplyr::select(subset, target_label, rho_disc = rho)
val1_dir <- val1 %>%
  dplyr::select(subset, target_label, rho_val = rho) %>%
  inner_join(disc_ref, by = c("subset", "target_label")) %>%
  mutate(consistent = sign(rho_disc) == sign(rho_val))
val2_dir <- val2 %>%
  dplyr::select(subset, target_label, rho_val = rho) %>%
  inner_join(disc_ref, by = c("subset", "target_label")) %>%
  mutate(consistent = sign(rho_disc) == sign(rho_val))

cat(sprintf("  Discovery vs GSE235236: %d / %d (%.0f%%)\n",
            sum(val1_dir$consistent), nrow(val1_dir),
            100*mean(val1_dir$consistent)))
cat(sprintf("  Discovery vs GSE75214 : %d / %d (%.0f%%)\n",
            sum(val2_dir$consistent), nrow(val2_dir),
            100*mean(val2_dir$consistent)))
combined_dir <- bind_rows(val1_dir %>% mutate(study="GSE235236"),
                           val2_dir %>% mutate(study="GSE75214"))
cat(sprintf("  Combined (both validations): %d / %d (%.0f%%)\n",
            sum(combined_dir$consistent), nrow(combined_dir),
            100*mean(combined_dir$consistent)))

# Per-target consistency
per_tg_cons <- combined_dir %>%
  group_by(target_label) %>%
  summarise(consistent = sum(consistent), total = n(),
            pct = 100*mean(consistent), .groups = "drop") %>%
  arrange(desc(pct))
print(as.data.frame(per_tg_cons))
write.csv(per_tg_cons, "results/11_direction_consistency_per_target.csv", row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Outputs:\n")
cat("  results/08_3cohort_effect_sizes.csv — 84 effect sizes\n")
cat("  results/09_per_cohort_pooled.csv\n")
cat("  results/10_3cohort_meta_pooled.csv ★ Pooled per target × scope\n")
cat("  results/11_direction_consistency_per_target.csv\n")
cat("  figures/06_3cohort_meta_forest.png ★ updated Fig 5C\n")
cat("  figures/07_3cohort_pooled_sensitivity.png ★ updated Fig 5D\n")
