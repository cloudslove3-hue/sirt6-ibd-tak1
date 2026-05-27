# ==============================================================================
# External validation cohort — GSE235236 (n=56: 8 Ctrl, 24 CD, 24 UC, sigmoid colon)
# Stage 1 파이프라인 적용: SIRT6 Spearman + WGCNA + (간단) partial correlation
# 목표: TAK1 negative coupling 재현 확인
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/external_validation")

for (d in c("results", "figures", "rds")) dir.create(d, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(ggplot2); library(ggpubr); library(pheatmap)
  library(DESeq2); library(WGCNA); library(matrixStats)
})
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))
set.seed(42)

# ==============================================================================
# 1. Metadata from series matrix
# ==============================================================================
cat(">>> [1/7] Parse series matrix\n")
sm <- readLines(gzfile("data/GSE235236_series_matrix.txt.gz"))
sample_rows <- sm[grepl("^!Sample_", sm)]

parse_row <- function(prefix) {
  ln <- sample_rows[startsWith(sample_rows, prefix)]
  if (!length(ln)) return(NULL)
  toks <- strsplit(sub(paste0("^", prefix, "\\s*"), "", ln[1]), "\t")[[1]]
  gsub('^"|"$', '', toks)
}
parse_characteristics <- function(key) {
  rows <- sample_rows[startsWith(sample_rows, "!Sample_characteristics_ch1")]
  for (r in rows) {
    toks <- gsub('^"|"$', '', strsplit(sub("^!Sample_characteristics_ch1\\s*", "", r), "\t")[[1]])
    if (length(toks) > 1 && startsWith(toks[1], paste0(key, ":"))) {
      return(sub(paste0("^", key, ":\\s*"), "", toks))
    }
  }
  NULL
}

gsm <- parse_row("!Sample_geo_accession")
title <- parse_row("!Sample_title")
disease <- parse_characteristics("disease state")
if (is.null(disease)) disease <- parse_characteristics("disease")
if (is.null(disease)) disease <- parse_characteristics("group")
location <- parse_characteristics("tissue")
if (is.null(location)) location <- parse_characteristics("anatomical site")
gender <- parse_characteristics("Sex")
if (is.null(gender)) gender <- parse_characteristics("gender")
age <- parse_characteristics("age")

# Inspect all characteristics keys
char_rows <- sample_rows[startsWith(sample_rows, "!Sample_characteristics_ch1")]
keys <- unique(sapply(char_rows, function(r) {
  toks <- gsub('^"|"$', '', strsplit(sub("^!Sample_characteristics_ch1\\s*", "", r), "\t")[[1]])
  sub(":.*", "", toks[1])
}))
cat("Characteristics keys found:", paste(keys, collapse = " | "), "\n")

pdata <- data.frame(
  gsm_id = gsm, title = title,
  disease = disease, location = location,
  gender = gender, age = age,
  stringsAsFactors = FALSE
)
cat("Sample distribution:\n")
print(table(pdata$disease))
write.csv(pdata, "results/00_pdata.csv", row.names = FALSE)

# ==============================================================================
# 2. Load all RSEM count files
# ==============================================================================
cat(">>> [2/7] Load 56 RSEM gene results files\n")
rsem_files <- list.files("data", pattern = "_RSEM\\.genes\\.results\\.gz$",
                         full.names = TRUE)
cat(sprintf("   files: %d\n", length(rsem_files)))

# 첫 파일로 gene_id 골격 확보
first <- fread(rsem_files[1], data.table = FALSE)
gene_ids <- first$gene_id   # ENSG#.#_SYMBOL

# Build expected_count matrix
counts <- sapply(rsem_files, function(f) {
  d <- fread(f, data.table = FALSE)
  d$expected_count[match(gene_ids, d$gene_id)]
})
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"
rownames(counts) <- gene_ids

# Filename → GSM mapping
fnames <- basename(rsem_files)
gsm_from_file <- sub("_.*", "", fnames)
colnames(counts) <- gsm_from_file

# Align with pdata
common <- intersect(colnames(counts), pdata$gsm_id)
counts <- counts[, common]
pdata <- pdata[match(common, pdata$gsm_id), ]
cat(sprintf("   matched samples: %d\n", length(common)))

# Symbol mapping (gene_id에서 추출)
symbol <- sub("^ENSG[0-9]+\\.[0-9]+_", "", rownames(counts))
ens_only <- sub("\\..*", "", sub("_.*", "", rownames(counts)))
# duplicate symbol → variance largest 유지
v <- matrixStats::rowVars(counts)
ord <- order(v, decreasing = TRUE)
counts <- counts[ord, ]; symbol <- symbol[ord]
keep <- !duplicated(symbol) & nchar(symbol) > 0
counts <- counts[keep, ]; symbol <- symbol[keep]
rownames(counts) <- symbol
cat(sprintf("   final counts: %d genes × %d samples\n",
            nrow(counts), ncol(counts)))

# ==============================================================================
# 3. QC
# ==============================================================================
cat(">>> [3/7] QC\n")
keep_g <- rowSums(counts >= 10) >= 0.1 * ncol(counts)
counts <- counts[keep_g, ]
cat(sprintf("   gene filter: %d\n", nrow(counts)))

lib <- colSums(counts)
keep_s <- lib >= 0.2 * median(lib)
counts <- counts[, keep_s]; pdata <- pdata[keep_s, ]
cat(sprintf("   QC pass: %d samples\n", ncol(counts)))

# disease 표준화
pdata$disease_clean <- case_when(
  grepl("^HC$|control|healthy|non.?IBD", pdata$disease, ignore.case = TRUE) ~ "Control",
  grepl("^UC$|ulcerative", pdata$disease, ignore.case = TRUE) ~ "UC",
  grepl("^CD$|crohn", pdata$disease, ignore.case = TRUE) ~ "CD",
  TRUE ~ pdata$disease
)
stopifnot(!any(is.na(pdata$disease_clean)))
pdata$disease_clean <- factor(pdata$disease_clean, levels = c("Control", "UC", "CD"))
cat("Disease (cleaned):\n"); print(table(pdata$disease_clean))

# ==============================================================================
# 4. DESeq2 vst
# ==============================================================================
cat(">>> [4/7] DESeq2 vst\n")
dds <- DESeqDataSetFromMatrix(counts, colData = pdata, design = ~ disease_clean)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)
saveRDS(expr, "rds/expr_vst.rds")
saveRDS(pdata, "rds/pdata_final.rds")
stopifnot("SIRT6" %in% rownames(expr))

# ==============================================================================
# 5. SIRT6 Spearman correlation
# ==============================================================================
cat(">>> [5/7] SIRT6 anchor Spearman\n")
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
  e <- expr[, samps]
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
write.csv(cor_res, "results/01_sirt6_correlation.csv", row.names = FALSE)
cat("Correlation table:\n"); print(as.data.frame(cor_res))

# Heatmap
cor_wide <- cor_res %>% dplyr::select(subset, target_label, rho) %>%
  pivot_wider(names_from = target_label, values_from = rho) %>%
  column_to_rownames("subset") %>% as.matrix()
pheatmap(cor_wide,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks = seq(-0.7, 0.7, length.out = 101),
  display_numbers = TRUE, fontsize_number = 10,
  main = "GSE235236 external — SIRT6 vs targets (Spearman ρ)",
  cluster_rows = FALSE, cluster_cols = TRUE,
  filename = "figures/01_validation_heatmap.png",
  width = 8, height = 4)

# ==============================================================================
# 6. WGCNA on Control + UC (n ~32)
# ==============================================================================
cat(">>> [6/7] WGCNA (focus subset)\n")
focus <- pdata$gsm_id[pdata$disease_clean %in% c("Control", "UC")]
v_genes <- order(matrixStats::rowVars(expr[, focus]), decreasing = TRUE)[1:5000]
sirt6_idx <- which(rownames(expr) == "SIRT6")
v_genes <- unique(c(v_genes, sirt6_idx))
expr_w <- t(expr[v_genes, focus])
cat(sprintf("   WGCNA input: %d samples × %d genes\n", nrow(expr_w), ncol(expr_w)))

if (nrow(expr_w) >= 15) {
  powers <- c(1:10, seq(12, 20, 2))
  sft <- pickSoftThreshold(expr_w, powerVector = powers, networkType = "signed",
                           verbose = 0)
  chosen_power <- if (is.na(sft$powerEstimate)) 12 else sft$powerEstimate
  cat(sprintf("   power = %d\n", chosen_power))

  net <- blockwiseModules(expr_w, power = chosen_power,
                          networkType = "signed", TOMType = "signed",
                          minModuleSize = 20, maxBlockSize = 6000,
                          mergeCutHeight = 0.25, numericLabels = TRUE,
                          verbose = 0, saveTOMs = FALSE)
  module_colors <- labels2colors(net$colors)
  names(module_colors) <- colnames(expr_w)
  sirt6_color <- unname(module_colors["SIRT6"])
  sirt6_module_label <- net$colors["SIRT6"]
  mod_genes <- names(net$colors)[net$colors == sirt6_module_label]
  cat(sprintf("   SIRT6 module color: %s, size: %d\n",
              sirt6_color, length(mod_genes)))

  ME <- net$MEs[, paste0("ME", sirt6_module_label)]
  kME <- apply(expr_w[, mod_genes], 2, function(x) cor(x, ME))
  hub <- data.frame(gene = names(kME), kME = kME) %>%
    arrange(desc(abs(kME))) %>% head(30)
  write.csv(hub, "results/02_sirt6_module_hub_top30.csv", row.names = FALSE)
  cat("   Top 30 hub:\n"); print(as.data.frame(hub))

  # Discovery hub와의 overlap
  disc_hub_file <- "../pilot_analysis/results/02_sirt6_module_hub_top30.csv"
  if (file.exists(disc_hub_file)) {
    disc_hub <- read.csv(disc_hub_file)
    overlap <- intersect(hub$gene, disc_hub$gene)
    cat(sprintf("\n   Hub gene overlap with GSE193677 (discovery): %d / 30\n",
                length(overlap)))
    cat("   shared:", paste(overlap, collapse = ", "), "\n")
    write.csv(data.frame(shared_hub = overlap),
              "results/03_hub_overlap_with_discovery.csv", row.names = FALSE)
  }

  saveRDS(net, "rds/wgcna_net.rds")
} else {
  cat("   Skip WGCNA: too few samples\n")
}

# ==============================================================================
# 7. Comparison with discovery cohort (GSE193677)
# ==============================================================================
cat(">>> [7/7] Compare with discovery cohort (GSE193677)\n")
disc_cor_file <- "../pilot_analysis/results/01_sirt6_correlation_table.csv"
if (file.exists(disc_cor_file)) {
  disc <- read.csv(disc_cor_file)
  # Match by subset (All_samples, Control, UC_all~UC, CD_all~CD)
  disc_match <- disc %>%
    mutate(subset = case_when(
      subset == "UC_all" ~ "UC",
      subset == "CD_all" ~ "CD",
      TRUE ~ subset
    )) %>%
    filter(subset %in% c("All_samples", "Control", "UC", "CD")) %>%
    dplyr::select(subset, target_label, rho_discovery = rho, q_discovery = q)
  cmp <- cor_res %>%
    dplyr::select(subset, target_label, rho_validation = rho, q_validation = q) %>%
    inner_join(disc_match, by = c("subset", "target_label")) %>%
    mutate(direction_consistent = sign(rho_discovery) == sign(rho_validation))
  print(as.data.frame(cmp))
  write.csv(cmp, "results/04_discovery_vs_validation.csv", row.names = FALSE)
  cat(sprintf("\n   Direction consistent: %d / %d\n",
              sum(cmp$direction_consistent, na.rm = TRUE), nrow(cmp)))

  # Scatter plot — discovery vs validation rho
  p_cmp <- ggplot(cmp, aes(rho_discovery, rho_validation,
                            color = subset, shape = target_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", alpha = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = target_label), vjust = -1, size = 3) +
    labs(title = "Discovery vs Validation — SIRT6 Spearman ρ",
         x = "GSE193677 (discovery, n=2,483)",
         y = "GSE235236 (validation, n≈56)") +
    theme_bw(base_size = 11)
  ggsave("figures/02_discovery_vs_validation.png", p_cmp,
         width = 8, height = 6, dpi = 130)
}

cat("\n=== DONE ===\n")
cat("Key result: TAK1(MAP3K7) Spearman ρ in validation =\n")
tak1_val <- cor_res %>% filter(target_label == "TAK1")
print(tak1_val)
