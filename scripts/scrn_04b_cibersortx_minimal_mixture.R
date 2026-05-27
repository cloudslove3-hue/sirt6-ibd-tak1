# ==============================================================================
# CIBERSORTx 업로드용 minimal bulk mixture
# 1.4GB 전체 파일은 압축 실패 / 업로드 한계 위험
# → signature matrix gene과 SC reference gene 합집합으로 subset (~200 gene)
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi); library(matrixStats)
})

# Signature matrix gene list
sig <- fread("results/07_CIBERSORTx_signature_matrix.tsv", data.table = FALSE)
sig_genes <- sig$Gene
cat("Signature genes:", length(sig_genes), "\n")

# Load full bulk count matrix
cat("Loading GSE193677 raw counts...\n")
counts_file <- "../pilot_analysis/data/GSE193677_counts.txt.gz"
hdr <- readLines(gzfile(counts_file), n = 1)
sample_ids <- gsub('^"|"$', '', strsplit(hdr, " ")[[1]])
counts_dt <- fread(counts_file, sep = " ", header = FALSE, skip = 1,
                   data.table = FALSE, quote = '"')
gene_ids <- counts_dt[, 1]
counts_bulk <- as.matrix(counts_dt[, -1])
storage.mode(counts_bulk) <- "integer"
rownames(counts_bulk) <- gene_ids
colnames(counts_bulk) <- sample_ids

# Ensembl → Symbol
ens_clean <- sub("\\..*$", "", rownames(counts_bulk))
sym_map <- mapIds(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
new_rn <- ifelse(is.na(sym_map), ens_clean, sym_map)
v <- matrixStats::rowVars(counts_bulk)
ord <- order(v, decreasing = TRUE)
counts_bulk <- counts_bulk[ord, ]; new_rn <- new_rn[ord]
counts_bulk <- counts_bulk[!duplicated(new_rn), ]
rownames(counts_bulk) <- new_rn[!duplicated(new_rn)]

# Subset to signature genes
common <- intersect(sig_genes, rownames(counts_bulk))
cat("Signature genes in bulk:", length(common), "\n")
counts_sub <- counts_bulk[common, ]

# CPM
cpm_sub <- edgeR::cpm(counts_sub, log = FALSE)
cat(sprintf("Subset mixture: %d genes × %d samples\n",
            nrow(cpm_sub), ncol(cpm_sub)))

bulk_df <- data.frame(Gene = rownames(cpm_sub), cpm_sub, check.names = FALSE)
write.table(bulk_df, "results/09b_CIBERSORTx_bulk_mixture_minimal.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
sz <- file.size("results/09b_CIBERSORTx_bulk_mixture_minimal.tsv") / 1024 / 1024
cat(sprintf("Saved: results/09b_CIBERSORTx_bulk_mixture_minimal.tsv (%.1f MB)\n", sz))

# 추가: full bulk를 청크 단위로 gzip 시도 (안전 fallback)
cat("\nAttempting safe gzip of full bulk mixture...\n")
full_file <- "results/09_CIBERSORTx_bulk_mixture.tsv"
gz_file <- "results/09_CIBERSORTx_bulk_mixture.tsv.gz"
if (file.exists(full_file)) {
  R.utils::gzip(full_file, destname = gz_file, overwrite = TRUE,
                remove = FALSE)
  sz2 <- file.size(gz_file) / 1024 / 1024
  cat(sprintf("Gzipped full bulk: %.1f MB\n", sz2))
}

cat("\n=== DONE ===\n")
cat("Recommended for CIBERSORTx upload:\n")
cat("  Signature: results/07_CIBERSORTx_signature_matrix.tsv (8.6 KB)\n")
cat("  Mixture: results/09b_CIBERSORTx_bulk_mixture_minimal.tsv (subset, signature genes only)\n")
cat("  또는 전체 mixture: results/09_CIBERSORTx_bulk_mixture.tsv.gz (gzipped)\n")
