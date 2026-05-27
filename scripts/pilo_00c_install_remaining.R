Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Library paths:\n"); print(.libPaths())

# 1. WGCNA — pure-R 패키지로 큰 어려움 없을 듯
cat("\n>>> Installing WGCNA...\n")
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA", lib = "C:/Rlibs", dependencies = TRUE)
}

# 2. DOSE 먼저 (clusterProfiler/enrichplot의 종속성)
cat("\n>>> Installing DOSE...\n")
if (!requireNamespace("DOSE", quietly = TRUE)) {
  BiocManager::install("DOSE", update = FALSE, ask = FALSE, lib = "C:/Rlibs")
}

# 3. enrichplot
cat("\n>>> Installing enrichplot...\n")
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot", update = FALSE, ask = FALSE, lib = "C:/Rlibs")
}

# 4. clusterProfiler
cat("\n>>> Installing clusterProfiler...\n")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE, lib = "C:/Rlibs")
}

# 5. matrixStats (used in script)
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  install.packages("matrixStats", lib = "C:/Rlibs")
}

# 검증
cat("\n=== FINAL CHECK ===\n")
miss <- c()
for (p in c("WGCNA", "DOSE", "enrichplot", "clusterProfiler", "matrixStats")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-20s : %s\n", p, if (ok) "OK" else "FAIL"))
  if (!ok) miss <- c(miss, p)
}
if (length(miss)) {
  cat("STILL MISSING:", paste(miss, collapse=", "), "\n")
  quit(status = 1)
} else {
  cat("\nAll required packages installed.\n")
}
