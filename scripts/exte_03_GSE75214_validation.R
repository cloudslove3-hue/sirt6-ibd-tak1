# ==============================================================================
# External validation cohort #2 — GSE75214 (Vancamelbeke et al. 2017)
# Platform: Affymetrix Human Gene 1.0 ST Array (GPL6244), microarray
# n = 194 (UC active 94 + inactive 16; CD colonic 8 + ileal active 51 + ileal inactive 16; control colon 11 + ileum 11)
# Tissue: colon + (neo-)terminal ileum
# Data: log2-normalized expression in series matrix
# Pipeline: Same Stage 1 (Spearman + WGCNA), platform-appropriate (no DESeq2 vst)
# 사전 등록: S1 Amendment 02 (2026-05-27)
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/external_validation")
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Affymetrix Hugene 1.0 ST annotation
if (!requireNamespace("hugene10sttranscriptcluster.db", quietly = TRUE))
  BiocManager::install("hugene10sttranscriptcluster.db",
                       update = FALSE, ask = FALSE, lib = "C:/Rlibs")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(ggplot2); library(ggpubr); library(pheatmap)
  library(WGCNA); library(matrixStats)
  library(hugene10sttranscriptcluster.db); library(AnnotationDbi)
})
set.seed(42)
for (d in c("results_GSE75214", "figures_GSE75214", "rds_GSE75214"))
  dir.create(d, showWarnings = FALSE)

# ==============================================================================
# 1. Parse series matrix
# ==============================================================================
cat(">>> [1/7] Parse GSE75214 series matrix\n")
sm_file <- "data_GSE75214/GSE75214_series_matrix.txt.gz"
sm <- readLines(gzfile(sm_file))
sample_rows <- sm[grepl("^!Sample_", sm)]

parse_row <- function(prefix) {
  ln <- sample_rows[startsWith(sample_rows, prefix)]
  if (!length(ln)) return(NULL)
  toks <- strsplit(sub(paste0("^", prefix, "\\s*"), "", ln[1]), "\t")[[1]]
  gsub('^"|"$', '', toks)
}
parse_chars <- function(key) {
  rows <- sample_rows[startsWith(sample_rows, "!Sample_characteristics_ch1")]
  for (r in rows) {
    toks <- gsub('^"|"$', '', strsplit(sub("^!Sample_characteristics_ch1\\s*", "", r), "\t")[[1]])
    if (length(toks) > 1 && startsWith(toks[1], paste0(key, ":")))
      return(sub(paste0("^", key, ":\\s*"), "", toks))
  }
  NULL
}

gsm <- parse_row("!Sample_geo_accession")
title <- parse_row("!Sample_title")
tissue <- parse_chars("tissue")
disease <- parse_chars("disease")
activity <- parse_chars("disease activity")
if (is.null(activity)) activity <- parse_chars("activity")

pdata <- data.frame(
  gsm_id = gsm, title = title, tissue = tissue,
  disease = disease, activity = activity,
  stringsAsFactors = FALSE
)
cat(sprintf("   n samples: %d\n", nrow(pdata)))
cat("   tissue × disease × activity:\n")
print(table(pdata$tissue, pdata$disease, pdata$activity, useNA = "ifany"))

# Disease cleaning
pdata$disease_clean <- case_when(
  grepl("control|healthy|normal", pdata$disease, ignore.case = TRUE) ~ "Control",
  grepl("ulcerative|^UC$", pdata$disease, ignore.case = TRUE) ~ "UC",
  grepl("crohn|^CD$", pdata$disease, ignore.case = TRUE) ~ "CD",
  TRUE ~ pdata$disease
)
pdata$disease_clean <- factor(pdata$disease_clean, levels = c("Control", "UC", "CD"))
cat("\nDisease (cleaned):\n"); print(table(pdata$disease_clean))

write.csv(pdata, "results_GSE75214/00_pdata.csv", row.names = FALSE)

# ==============================================================================
# 2. Load expression matrix
# ==============================================================================
cat(">>> [2/7] Load expression matrix from series matrix\n")
# Find table block start
all_lines <- readLines(gzfile(sm_file))
begin_idx <- which(all_lines == "!series_matrix_table_begin")
end_idx <- which(all_lines == "!series_matrix_table_end")
if (length(begin_idx) == 0) stop("table block not found")
expr_lines <- all_lines[(begin_idx + 1):(end_idx - 1)]
# write to temp + fread
tmp_file <- tempfile(fileext = ".txt")
writeLines(expr_lines, tmp_file)
expr_dt <- fread(tmp_file, sep = "\t", data.table = FALSE)
file.remove(tmp_file)

# First column = probe ID
probe_ids <- expr_dt[, 1]
expr_mat <- as.matrix(expr_dt[, -1])
rownames(expr_mat) <- as.character(probe_ids)
storage.mode(expr_mat) <- "numeric"
cat(sprintf("   expression: %d probes × %d samples\n", nrow(expr_mat), ncol(expr_mat)))

# Align with pdata
common <- intersect(colnames(expr_mat), pdata$gsm_id)
expr_mat <- expr_mat[, common]; pdata <- pdata[match(common, pdata$gsm_id), ]

# ==============================================================================
# 3. Probe → Gene symbol mapping
# ==============================================================================
cat(">>> [3/7] Affymetrix probe → SYMBOL mapping\n")
probe_to_sym <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = rownames(expr_mat),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
probe_to_sym <- probe_to_sym[!is.na(probe_to_sym$SYMBOL), ]
probe_to_sym <- probe_to_sym[!duplicated(probe_to_sym$PROBEID), ]
cat(sprintf("   probes with SYMBOL: %d / %d\n",
            nrow(probe_to_sym), nrow(expr_mat)))

# Map probe → symbol
sym_vec <- probe_to_sym$SYMBOL[match(rownames(expr_mat), probe_to_sym$PROBEID)]
keep <- !is.na(sym_vec)
expr_mat <- expr_mat[keep, ]; sym_vec <- sym_vec[keep]

# Multiple probes per gene → keep max-variance probe
v <- matrixStats::rowVars(expr_mat)
ord <- order(v, decreasing = TRUE)
expr_mat <- expr_mat[ord, ]; sym_vec <- sym_vec[ord]
keep_dedup <- !duplicated(sym_vec)
expr_mat <- expr_mat[keep_dedup, ]
rownames(expr_mat) <- sym_vec[keep_dedup]
cat(sprintf("   final symbol-level matrix: %d genes × %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

saveRDS(expr_mat, "rds_GSE75214/expr_GSE75214.rds")
saveRDS(pdata, "rds_GSE75214/pdata_GSE75214.rds")
stopifnot("SIRT6" %in% rownames(expr_mat))

# ==============================================================================
# 4. SIRT6 Spearman correlation per subgroup
# ==============================================================================
cat(">>> [4/7] SIRT6 anchor Spearman per subgroup\n")
targets <- c(NLRP3="NLRP3", CASP1="CASP1", IL1B="IL1B", ASC="PYCARD",
             TAK1="MAP3K7", NFKB_p65="RELA", FOXC1="FOXC1")

subgroups <- list(
  All_samples = pdata$gsm_id,
  Control     = pdata$gsm_id[pdata$disease_clean == "Control"],
  UC          = pdata$gsm_id[pdata$disease_clean == "UC"],
  CD          = pdata$gsm_id[pdata$disease_clean == "CD"]
)
cat("Subgroup sizes:\n")
for (n in names(subgroups)) cat(sprintf("   %s: n=%d\n", n, length(subgroups[[n]])))

compute_cor <- function(samps, label) {
  if (length(samps) < 5) return(NULL)
  e <- expr_mat[, samps]
  s <- e["SIRT6", ]
  bind_rows(lapply(names(targets), function(tg) {
    g <- targets[tg]
    if (!g %in% rownames(e)) return(NULL)
    ct <- cor.test(s, e[g, ], method = "spearman", exact = FALSE)
    data.frame(subset = label, n = length(samps),
               target_label = tg, target_gene = unname(g),
               rho = unname(ct$estimate), p = ct$p.value)
  }))
}
cor_res <- bind_rows(lapply(names(subgroups), function(n)
  compute_cor(subgroups[[n]], n))) %>%
  group_by(subset) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
write.csv(cor_res, "results_GSE75214/01_sirt6_correlation.csv", row.names = FALSE)
cat("Correlation table:\n"); print(as.data.frame(cor_res))

# Heatmap
cor_wide <- cor_res %>% dplyr::select(subset, target_label, rho) %>%
  pivot_wider(names_from = target_label, values_from = rho) %>%
  column_to_rownames("subset") %>% as.matrix()
pheatmap(cor_wide,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks = seq(-0.7, 0.7, length.out = 101),
  display_numbers = TRUE, fontsize_number = 10,
  main = "GSE75214 (n=194, Vancamelbeke 2017) — SIRT6 vs targets",
  cluster_rows = FALSE, cluster_cols = TRUE,
  filename = "figures_GSE75214/01_GSE75214_heatmap.png",
  width = 8, height = 4)

# ==============================================================================
# 5. WGCNA (n=194 충분)
# ==============================================================================
cat(">>> [5/7] WGCNA on UC+Control subset\n")
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))
# AnnotationDbi::select 등이 cor()을 마스킹할 수 있으므로 명시
cor_orig <- cor
cor <- WGCNA::cor
on.exit(assign("cor", cor_orig, envir = .GlobalEnv), add = TRUE)
focus <- pdata$gsm_id[pdata$disease_clean %in% c("Control", "UC")]
cat(sprintf("   focus subset: %d samples\n", length(focus)))

v_genes <- order(matrixStats::rowVars(expr_mat[, focus]), decreasing = TRUE)[1:5000]
sirt6_idx <- which(rownames(expr_mat) == "SIRT6")
v_genes <- unique(c(v_genes, sirt6_idx))
expr_w <- t(expr_mat[v_genes, focus])
cat(sprintf("   WGCNA input: %d × %d\n", nrow(expr_w), ncol(expr_w)))

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(expr_w, powerVector = powers, networkType = "signed",
                         verbose = 0)
chosen_power <- if (is.na(sft$powerEstimate)) 12 else sft$powerEstimate
cat(sprintf("   power = %d\n", chosen_power))

net <- blockwiseModules(expr_w, power = chosen_power,
                        networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, maxBlockSize = 6000,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        verbose = 0, saveTOMs = FALSE)
module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_w)
sirt6_color <- unname(module_colors["SIRT6"])
sirt6_module_label <- net$colors["SIRT6"]
mod_genes <- names(net$colors)[net$colors == sirt6_module_label]
cat(sprintf("   SIRT6 module: %s, size: %d\n", sirt6_color, length(mod_genes)))

ME <- net$MEs[, paste0("ME", sirt6_module_label)]
kME <- apply(expr_w[, mod_genes], 2, function(x) cor(x, ME))
hub <- data.frame(gene = names(kME), kME = kME) %>%
  arrange(desc(abs(kME))) %>% head(30)
write.csv(hub, "results_GSE75214/02_sirt6_module_hub_top30.csv", row.names = FALSE)
cat("Top 30 hub:\n"); print(as.data.frame(hub))

# Compare with discovery hub
disc_hub_file <- "../pilot_analysis/results/02_sirt6_module_hub_top30.csv"
if (file.exists(disc_hub_file)) {
  disc_hub <- read.csv(disc_hub_file)
  overlap <- intersect(hub$gene, disc_hub$gene)
  cat(sprintf("\n   Hub overlap with discovery: %d / 30\n", length(overlap)))
  cat("   shared:", paste(overlap, collapse = ", "), "\n")
}
saveRDS(net, "rds_GSE75214/wgcna_net.rds")

# ==============================================================================
# 6. Comparison with discovery cohort
# ==============================================================================
cat(">>> [6/7] Compare with discovery (GSE193677)\n")
disc_cor_file <- "../pilot_analysis/results/01_sirt6_correlation_table.csv"
if (file.exists(disc_cor_file)) {
  disc <- read.csv(disc_cor_file) %>%
    mutate(subset = case_when(
      subset == "UC_all" ~ "UC", subset == "CD_all" ~ "CD",
      TRUE ~ subset)) %>%
    filter(subset %in% c("All_samples", "Control", "UC", "CD")) %>%
    dplyr::select(subset, target_label, rho_disc = rho, q_disc = q)
  cmp <- cor_res %>%
    dplyr::select(subset, target_label, rho_val = rho, q_val = q) %>%
    inner_join(disc, by = c("subset", "target_label")) %>%
    mutate(direction_consistent = sign(rho_disc) == sign(rho_val))
  write.csv(cmp, "results_GSE75214/03_disc_vs_GSE75214.csv", row.names = FALSE)
  print(as.data.frame(cmp))
  cat(sprintf("\n   Direction consistent: %d / %d (%.0f%%)\n",
              sum(cmp$direction_consistent), nrow(cmp),
              100*mean(cmp$direction_consistent)))
}

cat("\n=== DONE ===\n")
cat("Key results:\n")
tak1 <- cor_res %>% filter(target_label == "TAK1")
print(as.data.frame(tak1))
