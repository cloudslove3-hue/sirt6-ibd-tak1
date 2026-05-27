# ==============================================================================
# 사전 단계: 모든 필수 패키지 설치
# 환경: R 4.6.0 / Windows
# ==============================================================================

# 0. 개인 라이브러리 설정 (관리자 권한 없이 설치)
user_lib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library", "4.6")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
cat("Library paths:\n"); print(.libPaths())

# CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# 1. BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}
BiocManager::install(update = FALSE, ask = FALSE, lib = user_lib)

# 2. CRAN 패키지 일괄
cran_pkgs <- c(
  "tidyverse", "data.table", "ggplot2", "ggpubr", "corrplot",
  "pheatmap", "RColorBrewer", "ggrepel", "rmarkdown", "knitr",
  "DT", "R.utils"
)
cat(">>> Installing CRAN packages...\n")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("  Installing", p, "...\n")
    install.packages(p, lib = user_lib, dependencies = TRUE)
  } else {
    cat("  Already installed:", p, "\n")
  }
}

# 3. Bioconductor 패키지
bioc_pkgs <- c(
  "GEOquery", "DESeq2", "limma", "edgeR", "sva",
  "WGCNA", "clusterProfiler", "org.Hs.eg.db",
  "AnnotationDbi", "ReactomePA"
)
cat(">>> Installing Bioconductor packages...\n")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("  Installing", p, "...\n")
    BiocManager::install(p, update = FALSE, ask = FALSE, lib = user_lib)
  } else {
    cat("  Already installed:", p, "\n")
  }
}

# 4. 설치 검증
cat("\n=== Verification ===\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing <- c()
for (p in all_pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-25s : %s\n", p, if (ok) "OK" else "FAIL"))
  if (!ok) missing <- c(missing, p)
}

if (length(missing) > 0) {
  cat("\nMISSING:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("\nAll packages installed successfully.\n")
}
