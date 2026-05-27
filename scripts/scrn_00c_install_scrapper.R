Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# scrapper — new SingleR backend
if (!requireNamespace("scrapper", quietly = TRUE)) {
  cat("Installing scrapper...\n")
  tryCatch(
    install.packages("scrapper", lib = "C:/Rlibs", dependencies = TRUE),
    error = function(e) {
      cat("CRAN failed, trying Bioconductor...\n")
      BiocManager::install("scrapper", update = FALSE, ask = FALSE, lib = "C:/Rlibs")
    }
  )
}

# alabaster.base could also be needed
if (!requireNamespace("alabaster.base", quietly = TRUE)) {
  BiocManager::install("alabaster.base", update = FALSE, ask = FALSE, lib = "C:/Rlibs")
}

for (p in c("scrapper", "alabaster.base", "SingleR")) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-15s : %s\n", p, if (ok) "OK" else "FAIL"))
}
