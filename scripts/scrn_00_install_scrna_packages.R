# ==============================================================================
# scRNA 분석을 위한 추가 패키지 설치
# 기반: 기존 pilot_analysis 환경 (R 4.6.0, C:/Rlibs)
# ==============================================================================
Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c("Seurat", "SeuratObject", "patchwork", "remotes",
               "harmony", "scCustomize", "viridis", "cowplot", "ggrepel",
               "BisqueRNA")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing CRAN:", p, "\n")
    install.packages(p, lib = "C:/Rlibs", dependencies = TRUE)
  }
}

bioc_pkgs <- c("SingleCellExperiment", "scran", "scater", "scuttle",
               "Biobase", "celldex", "SingleR", "DropletUtils",
               "MatrixGenerics", "S4Vectors")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing Bioc:", p, "\n")
    tryCatch(
      BiocManager::install(p, update = FALSE, ask = FALSE, lib = "C:/Rlibs"),
      error = function(e) cat("  ERROR:", p, ":", conditionMessage(e), "\n")
    )
  }
}

# MuSiC from GitHub
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  cat("Installing MuSiC from GitHub...\n")
  remotes::install_github("xuranw/MuSiC", upgrade = "never", lib = "C:/Rlibs")
}

cat("\n=== Verification ===\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs, "MuSiC")
miss <- c()
for (p in all_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-22s : %s\n", p, if (ok) "OK" else "FAIL"))
  if (!ok) miss <- c(miss, p)
}
if (length(miss)) {
  cat("MISSING:", paste(miss, collapse=", "), "\n")
  quit(status = 1)
}
cat("All packages OK.\n")
