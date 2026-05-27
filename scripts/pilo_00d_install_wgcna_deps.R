Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Installing WGCNA Bioconductor dependencies...\n")
for (p in c("impute", "preprocessCore", "minet", "GO.db")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(">>>", p, "\n")
    BiocManager::install(p, update = FALSE, ask = FALSE, lib = "C:/Rlibs")
  } else {
    cat("OK:", p, "\n")
  }
}

# minet은 archived(retired)일 수 있음 — 없으면 무시
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  cat("Re-installing WGCNA after dependencies...\n")
  install.packages("WGCNA", lib = "C:/Rlibs", dependencies = FALSE)
}

cat("\n=== CHECK ===\n")
for (p in c("impute", "preprocessCore", "minet", "GO.db", "WGCNA",
            "DOSE", "enrichplot", "clusterProfiler", "matrixStats",
            "GEOquery", "DESeq2", "limma", "edgeR", "sva",
            "org.Hs.eg.db", "AnnotationDbi", "data.table", "tidyverse")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-22s : %s\n", p, if (ok) "OK" else "FAIL"))
}
