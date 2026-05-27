Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# metafor — random-effects meta-analysis
# ppcor — partial correlation sensitivity check (우선순위 4d)
pkgs <- c("metafor", "ppcor")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing", p, "\n")
    install.packages(p, lib = "C:/Rlibs", dependencies = TRUE)
  }
}
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-10s : %s\n", p, if (ok) "OK" else "FAIL"))
}
