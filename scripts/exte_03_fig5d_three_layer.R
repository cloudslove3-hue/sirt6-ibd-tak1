# ==============================================================================
# Fig 5D v2 — 3-layer pooled effect sizes
# Discovery-only / Validation-only / All-pooled per target
# Reviewer "discovery weight dominance" 사전 방어
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/external_validation")

suppressPackageStartupMessages({
  library(metafor); library(dplyr); library(tidyr); library(ggplot2)
  library(tibble); library(stringr)
})

# Load effect sizes
es <- read.csv("results/05_meta_effect_sizes.csv")
cat(sprintf("Loaded %d effect sizes\n", nrow(es)))

# Per target × scope (Discovery / Validation / All) → pooled REML
targets <- unique(es$target_label)
scopes <- list(
  Discovery   = "Discovery (GSE193677)",
  Validation  = "Validation (GSE235236)",
  All_studies = c("Discovery (GSE193677)", "Validation (GSE235236)")
)

meta_3layer <- list()
for (tg in targets) {
  for (sc in names(scopes)) {
    sub <- es %>% filter(target_label == tg, study %in% scopes[[sc]])
    if (nrow(sub) < 2) next
    m <- rma(yi = sub$yi, vi = sub$vi, method = "REML")
    z <- as.numeric(m$beta)
    meta_3layer[[length(meta_3layer)+1]] <- data.frame(
      target = tg, scope = sc, k = nrow(sub),
      pooled_z = z,
      pooled_rho   = (exp(2*z)-1)/(exp(2*z)+1),
      pooled_lower = (exp(2*m$ci.lb)-1)/(exp(2*m$ci.lb)+1),
      pooled_upper = (exp(2*m$ci.ub)-1)/(exp(2*m$ci.ub)+1),
      p_pooled = m$pval,
      I2 = m$I2,
      tau2 = m$tau2
    )
  }
}
meta_df <- bind_rows(meta_3layer)
write.csv(meta_df, "results/07_meta_3layer_pooled.csv", row.names = FALSE)
print(meta_df)

# Order targets by All-studies pooled ρ
target_order <- meta_df %>% filter(scope == "All_studies") %>%
  arrange(pooled_rho) %>% pull(target)
meta_df$target <- factor(meta_df$target, levels = rev(target_order))
meta_df$scope  <- factor(meta_df$scope,
                         levels = c("Discovery", "Validation", "All_studies"))

# Y position offset within each target for the 3 layers
meta_df$y_offset <- as.numeric(meta_df$scope) * 0.22 - 0.44
meta_df$y_pos    <- as.numeric(meta_df$target) + meta_df$y_offset

# I² label only for All_studies row
meta_df$label_I2 <- ifelse(meta_df$scope == "All_studies",
                           sprintf("I²=%.0f%%, k=%d", meta_df$I2, meta_df$k),
                           "")

p <- ggplot(meta_df, aes(y = y_pos, x = pooled_rho, color = scope)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = pooled_lower, xmax = pooled_upper),
                 height = 0.12, linewidth = 0.6) +
  geom_point(aes(shape = scope, size = scope)) +
  geom_text(aes(label = label_I2, x = pmax(pooled_upper + 0.05, 0.55)),
            hjust = 0, size = 3.1, color = "gray20") +
  scale_y_continuous(breaks = seq_along(levels(meta_df$target)),
                     labels = levels(meta_df$target),
                     expand = expansion(add = 0.5)) +
  scale_x_continuous(limits = c(-1, 1.2),
                     breaks = seq(-1, 1, 0.25)) +
  scale_color_manual(values = c(Discovery = "#1f78b4",
                                Validation = "#e31a1c",
                                All_studies = "black")) +
  scale_shape_manual(values = c(Discovery = 16, Validation = 17, All_studies = 18)) +
  scale_size_manual(values = c(Discovery = 3, Validation = 3, All_studies = 4.5)) +
  labs(x = "Pooled Spearman ρ (random-effects REML, 95% CI)",
       y = NULL, color = NULL, shape = NULL, size = NULL,
       title = "Cross-cohort pooled effect sizes — 3-layer",
       subtitle = "Discovery-only (blue), Validation-only (red), All pooled (black). I² shown for all-studies pooled.") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("figures/05_meta_pooled_3layer.png", p,
       width = 10.5, height = 5.5, dpi = 150)

cat("\n=== DONE ===\n")
cat("  results/07_meta_3layer_pooled.csv\n")
cat("  figures/05_meta_pooled_3layer.png ★ replaces Fig 5D\n")
