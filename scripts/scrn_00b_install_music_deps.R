Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# TOAST — Bioconductor
if (!requireNamespace("TOAST", quietly = TRUE)) {
  BiocManager::install("TOAST", update = FALSE, ask = FALSE, lib = "C:/Rlibs")
}

# MuSiC 재시도
if (!requireNamespace("MuSiC", quietly = TRUE)) {
  remotes::install_github("xuranw/MuSiC", upgrade = "never", lib = "C:/Rlibs",
                          dependencies = TRUE)
}

# BisqueRNA 재시도
if (!requireNamespace("BisqueRNA", quietly = TRUE)) {
  install.packages("BisqueRNA", lib = "C:/Rlibs", dependencies = TRUE)
}

for (p in c("TOAST", "MuSiC", "BisqueRNA")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-15s : %s\n", p, if (ok) "OK" else "FAIL"))
}
