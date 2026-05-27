# ==============================================================================
# Triple validation of partial Spearman
# (1) rank-residual lm (본 분석) — 이미 계산됨
# (2) ppcor::pcor.test (method="spearman") — 이미 계산됨, r=1.0000으로 일치
# (3) psych::partial.r on Spearman correlation matrix — 독립 제3 구현
# (4) 손계산: 정확한 공식 ρ_xy.z = (ρ_xy − ρ_xz·ρ_yz) / sqrt((1−ρ_xz²)(1−ρ_yz²))
#
# 목적: ppcor와의 r=1.0000이 algebraic equivalence임을 확인.
# psych::partial.r과 손계산은 inverse-correlation-matrix 공식에 기반 →
# rank-residual lm과 동일한 결과가 나와야 algebraic 일치 확정.
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

# Install psych if needed
if (!requireNamespace("psych", quietly = TRUE))
  install.packages("psych", lib = "C:/Rlibs", repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ppcor); library(psych); library(dplyr); library(tidyr)
})

# Load data
expr <- readRDS("../pilot_analysis/rds/expr_vst.rds")
pdata <- readRDS("../pilot_analysis/rds/pdata_final.rds")
prop_est <- read.csv("results/04_MuSiC_cell_fractions.csv")
common <- intersect(colnames(expr), prop_est$sample)
expr <- expr[, common]
pdata_sub <- pdata[match(common, pdata$biopsy_id), ]
epi_col <- grep("epi", colnames(prop_est), ignore.case = TRUE, value = TRUE)[1]
epi_frac <- prop_est[[epi_col]][match(common, prop_est$sample)]

# Test cases: spot check 4 representative subgroup×target combinations
cases <- list(
  list(name = "All_samples × TAK1 (MAP3K7)", samps = common, target = "MAP3K7"),
  list(name = "UC_inflamed × NLRP3",
       samps = common[pdata_sub$disease == "UC" & pdata_sub$inflammation == "I"],
       target = "NLRP3"),
  list(name = "CD_inflamed × RELA",
       samps = common[pdata_sub$disease == "CD" & pdata_sub$inflammation == "I"],
       target = "RELA"),
  list(name = "Control × FOXC1",
       samps = common[pdata_sub$disease == "Control"], target = "FOXC1")
)

cat("=== Triple-implementation validation ===\n\n")
results <- list()
for (cs in cases) {
  idx <- match(cs$samps, common)
  x <- expr["SIRT6", idx]
  y <- expr[cs$target, idx]
  z <- epi_frac[idx]

  # (1) rank-residual lm (본 분석)
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  r1 <- cor(resid(lm(rx ~ rz)), resid(lm(ry ~ rz)), method = "pearson")

  # (2) ppcor::pcor.test method="spearman"
  r2 <- ppcor::pcor.test(x, y, z, method = "spearman")$estimate

  # (3) psych::partial.r — ranks의 Spearman matrix → partial
  ranks_mat <- cbind(rx, ry, rz)
  rho_mat <- cor(ranks_mat, method = "pearson")  # ranks의 Pearson = Spearman
  r3 <- psych::partial.r(rho_mat, x = c("rx","ry"), y = "rz", method = "pearson")["rx","ry"]

  # (4) 손계산: ρ_xy.z formula
  rho_xy <- cor(rx, ry, method = "pearson")
  rho_xz <- cor(rx, rz, method = "pearson")
  rho_yz <- cor(ry, rz, method = "pearson")
  r4 <- (rho_xy - rho_xz * rho_yz) /
        sqrt((1 - rho_xz^2) * (1 - rho_yz^2))

  cat(sprintf("[%s]\n", cs$name))
  cat(sprintf("  (1) rank-residual lm:      %.10f\n", r1))
  cat(sprintf("  (2) ppcor::pcor.test:      %.10f\n", r2))
  cat(sprintf("  (3) psych::partial.r:      %.10f\n", r3))
  cat(sprintf("  (4) closed-form formula:   %.10f\n", r4))
  cat(sprintf("  max pairwise diff:         %.2e\n\n",
              max(abs(c(r1,r2,r3,r4) - r1))))

  results[[length(results)+1]] <- data.frame(
    case = cs$name,
    rankres = r1, ppcor = r2, psych = r3, formula = r4,
    max_diff = max(abs(c(r1,r2,r3,r4) - r1))
  )
}
res_df <- bind_rows(results)
write.csv(res_df, "results/15_partial_triple_validation.csv", row.names = FALSE)

cat("\n=== CONCLUSION ===\n")
if (max(res_df$max_diff) < 1e-8) {
  cat("All four implementations agree within floating-point precision.\n")
  cat("→ ALGEBRAICALLY EQUIVALENT (사용자 의심 확정)\n")
  cat("→ Methods 문구를 'algebraically equivalent for rank-transformed inputs and\n")
  cat("   yielded numerically indistinguishable estimates within floating-point precision'\n")
  cat("   으로 교체 권장. 별도 sensitivity check (Bootstrap CI 등)로 보강 필요.\n")
} else {
  cat("Implementations differ. Non-equivalent estimators.\n")
}
