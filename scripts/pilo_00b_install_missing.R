# 한글 사용자명 호환을 위해 TMP를 ASCII 경로로 강제
Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
dir.create("C:/Temp/R", showWarnings = FALSE, recursive = TRUE)

# 라이브러리도 ASCII 경로 사용
user_lib <- "C:/Rlibs"
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

cat("Library paths:\n"); print(.libPaths())
cat("Temp dir:", tempdir(), "\n")

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = user_lib)

# 재설치 필요 — 그리고 이전에 한글 경로에 설치된 패키지도 ASCII lib로 통합 설치
all_needed <- c(
  "tidyverse", "data.table", "ggplot2", "ggpubr", "corrplot",
  "pheatmap", "RColorBrewer", "ggrepel", "rmarkdown", "knitr",
  "DT", "R.utils",
  "GEOquery", "DESeq2", "limma", "edgeR", "sva",
  "WGCNA", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
)
# ReactomePA는 enrichplot 의존성 — Reactome 우회 위해 skip 가능. 일단 시도.

bioc_set <- c("GEOquery", "DESeq2", "limma", "edgeR", "sva",
              "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
cran_set <- setdiff(all_needed, bioc_set)

cat(">>> CRAN packages\n")
for (p in cran_set) {
  if (!requireNamespace(p, quietly = TRUE, lib.loc = user_lib)) {
    cat("  Installing", p, "...\n")
    tryCatch(
      install.packages(p, lib = user_lib, dependencies = TRUE),
      error = function(e) cat("  ERROR:", p, ":", conditionMessage(e), "\n")
    )
  } else {
    cat("  OK:", p, "\n")
  }
}

cat(">>> Bioconductor packages\n")
for (p in bioc_set) {
  if (!requireNamespace(p, quietly = TRUE, lib.loc = user_lib)) {
    cat("  Installing", p, "...\n")
    tryCatch(
      BiocManager::install(p, update = FALSE, ask = FALSE, lib = user_lib),
      error = function(e) cat("  ERROR:", p, ":", conditionMessage(e), "\n")
    )
  } else {
    cat("  OK:", p, "\n")
  }
}

cat("\n=== Final verification ===\n")
missing <- c()
for (p in all_needed) {
  ok <- requireNamespace(p, quietly = TRUE, lib.loc = user_lib)
  cat(sprintf("  %-20s : %s\n", p, if (ok) "OK" else "FAIL"))
  if (!ok) missing <- c(missing, p)
}
if (length(missing) > 0) {
  cat("\nSTILL MISSING:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("\nAll packages OK.\n")
}
