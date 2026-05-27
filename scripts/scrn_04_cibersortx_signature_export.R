# ==============================================================================
# CIBERSORTx-compatible signature matrix export (cross-validation 용)
# MuSiC을 main으로, CIBERSORTx로 cross-validate.
# Smillie 51 cell type → broad 4 category(Epi/Immune/Stromal/Other) 압축 후
# CIBERSORTx 웹 도구에 업로드할 수 있는 reference matrix 생성.
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(scuttle)
  library(dplyr); library(tidyr); library(tibble); library(Matrix)
})

# ==============================================================================
# 1. Load annotated Seurat
# ==============================================================================
cat(">>> [1/4] Load Smillie Seurat\n")
seu <- readRDS("rds/smillie_seurat_annotated.rds")
cat(sprintf("   cells: %d, genes: %d\n", ncol(seu), nrow(seu)))
cat("   broad_type:\n"); print(table(seu$broad_type))

# ==============================================================================
# 2. Marker gene 동정 (cell-type-specific)
# ==============================================================================
cat(">>> [2/4] Identify cell-type marker genes\n")
Idents(seu) <- "broad_type"
markers <- FindAllMarkers(seu,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.5,
                          test.use = "wilcox",
                          verbose = FALSE)
write.csv(markers, "results/06_celltype_markers.csv", row.names = FALSE)

# Top 50 markers per cell type
top50 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)
cat(sprintf("   %d cell types, %d total markers (top 50 each)\n",
            length(unique(top50$cluster)), nrow(top50)))

# ==============================================================================
# 3. CIBERSORTx-format reference signature matrix
# ==============================================================================
cat(">>> [3/4] Build CIBERSORTx-style signature matrix\n")
# Format: rows = genes, columns = cell types, values = mean log-normalized expression
sig_genes <- unique(top50$gene)

expr_norm <- GetAssayData(seu, assay = "RNA", layer = "data")
expr_sig <- expr_norm[intersect(sig_genes, rownames(expr_norm)), ]

# Mean expression per broad cell type
broad_factor <- factor(seu$broad_type)
sig_matrix <- sapply(levels(broad_factor), function(ct) {
  cells_ct <- which(broad_factor == ct)
  rowMeans(expr_sig[, cells_ct, drop = FALSE])
})

# CIBERSORTx 기본 input format: TSV, first column = gene symbol, no row.names index
sig_df <- as.data.frame(sig_matrix) %>%
  rownames_to_column("Gene")
write.table(sig_df, "results/07_CIBERSORTx_signature_matrix.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("   signature matrix: %d genes × %d cell types\n",
            nrow(sig_df), ncol(sig_df) - 1))

# ==============================================================================
# 4. Single-cell raw counts (CIBERSORTx Build Mode 옵션용)
# ==============================================================================
cat(">>> [4/4] Export single-cell reference (CIBERSORTx Build Mode)\n")
# CIBERSORTx Single-Cell Profiles Mode 입력:
#   - row 1: cell type labels
#   - row 2+: gene symbol + raw count
# Smillie cell 11K은 가능

counts_sc <- GetAssayData(seu, assay = "RNA", layer = "counts")
# Sparse → dense (CIBERSORTx 요구). 큰 문제 안 됨 — 11K cells.
counts_dense <- as.matrix(counts_sc)

# Header: gene + cell types
ct_labels <- seu$broad_type
sc_df <- data.frame(Gene = rownames(counts_dense), counts_dense, check.names = FALSE)
colnames(sc_df)[-1] <- ct_labels  # CIBERSORTx 요구: 열 헤더가 cell type label

write.table(sc_df, "results/08_CIBERSORTx_SC_reference.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("   single-cell reference: %d genes × %d cells\n",
            nrow(sc_df), ncol(sc_df) - 1))

# ==============================================================================
# 5. Bulk mixture export (CIBERSORTx 업로드용)
# ==============================================================================
# CIBERSORTx는 raw count가 아닌 TPM/CPM 권장. vst(log-transformed) 사용 불가.
# raw count → CPM 변환 후 export
cat(">>> [Bonus] Export GSE193677 as CIBERSORTx mixture file (CPM)\n")
library(data.table); library(edgeR)
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
library(org.Hs.eg.db); library(AnnotationDbi); library(matrixStats)
ens_clean <- sub("\\..*$", "", rownames(counts_bulk))
sym_map <- mapIds(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
new_rn <- ifelse(is.na(sym_map), ens_clean, sym_map)
v <- matrixStats::rowVars(counts_bulk)
ord <- order(v, decreasing = TRUE)
counts_bulk <- counts_bulk[ord, ]; new_rn <- new_rn[ord]
counts_bulk <- counts_bulk[!duplicated(new_rn), ]
rownames(counts_bulk) <- new_rn[!duplicated(new_rn)]

# CPM (counts per million)
cpm_bulk <- edgeR::cpm(counts_bulk, log = FALSE)
bulk_df <- data.frame(Gene = rownames(cpm_bulk), cpm_bulk, check.names = FALSE)
write.table(bulk_df, "results/09_CIBERSORTx_bulk_mixture.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("   bulk mixture (CPM): %d genes × %d samples\n",
            nrow(bulk_df), ncol(bulk_df) - 1))

cat("\n=== CIBERSORTx export DONE ===\n")
cat("CIBERSORTx 웹(https://cibersortx.stanford.edu/) 업로드:\n")
cat("  Step 1: Custom Module → 'Create Signature Matrix'\n")
cat("    - Single Cell File: results/08_CIBERSORTx_SC_reference.tsv\n")
cat("    또는 이미 signature 있으면:\n")
cat("    - Signature Matrix File: results/07_CIBERSORTx_signature_matrix.tsv\n")
cat("  Step 2: 'Impute Cell Fractions'\n")
cat("    - Signature Matrix: 위에서 생성한 것\n")
cat("    - Mixture File: results/09_CIBERSORTx_bulk_mixture.tsv\n")
cat("    - Batch correction (S-mode) ON\n")
cat("    - Permutations: 100\n")
cat("  결과를 MuSiC 결과와 cross-validate.\n")
