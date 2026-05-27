# ==============================================================================
# SIRT6 1차 시범 분석 v2 — GSE193677 (Mount Sinai IBD cohort)
# - 데이터/메타데이터 구조 실측 확인 후 작성된 ready-to-run 버전
# - 사전 다운로드된 파일 사용: data/GSE193677_counts.txt.gz, data/GSE193677_series_matrix.txt.gz
# - 한글 사용자명 호환을 위해 TMP/lib는 ASCII 경로
# ==============================================================================

# ---- 환경 설정 ----
Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
dir.create("C:/Temp/R", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("C:/Rlibs", .libPaths()))

setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/pilot_analysis")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("rds", showWarnings = FALSE)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(matrixStats)
  library(ggplot2); library(ggpubr); library(pheatmap)
  library(DESeq2); library(WGCNA); library(AnnotationDbi); library(org.Hs.eg.db)
  library(clusterProfiler)
})
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))

# ==============================================================================
# 1. 메타데이터 파싱 — series matrix 직접 처리 (GEOquery 우회)
# ==============================================================================
cat(">>> [1/8] Parsing series matrix...\n")
sm_file <- "data/GSE193677_series_matrix.txt.gz"
sm <- readLines(gzfile(sm_file))
# !Sample_* 행 추출
sample_rows <- sm[grepl("^!Sample_", sm)]

parse_row <- function(row_prefix) {
  line <- sample_rows[startsWith(sample_rows, row_prefix)]
  if (length(line) == 0) return(NULL)
  # 첫 번째만 사용 (characteristics는 여러 줄 있을 수 있음)
  line[1] %>%
    sub(paste0("^", row_prefix, "\\s*"), "", .) %>%
    strsplit("\t") %>% .[[1]] %>%
    str_replace_all('^"|"$', "")
}

parse_characteristics <- function(key) {
  # !Sample_characteristics_ch1 행 중 첫 토큰이 "key:"인 것 찾기
  rows <- sample_rows[startsWith(sample_rows, "!Sample_characteristics_ch1")]
  for (r in rows) {
    toks <- r %>% sub("^!Sample_characteristics_ch1\\s*", "", .) %>%
            strsplit("\t") %>% .[[1]] %>% str_replace_all('^"|"$', "")
    if (length(toks) > 1 && startsWith(toks[1], paste0(key, ":"))) {
      return(sub(paste0("^", key, ":\\s*"), "", toks))
    }
  }
  NULL
}

titles <- parse_row("!Sample_title")
gsm    <- parse_row("!Sample_geo_accession")

# 사용자 정의 characteristics
disease   <- parse_characteristics("ibd_disease")          # UC / CD / Control
location  <- parse_characteristics("regionre")             # Rectum, Ileum, ...
inflam    <- parse_characteristics("typere")               # I / NonI
disease2  <- parse_characteristics("diseasetypere")        # UC.I / UC.NonI / CD.I / ...
age       <- parse_characteristics("study_eligibility_age_at_endo")
gender    <- parse_characteristics("demographics_gender")
endo_sev  <- parse_characteristics("ibd_endoseverity_4levels")  # Inactive/Mild/Mod/Severe
nancy     <- parse_characteristics("nancyindex")
mayo_uc   <- parse_characteristics("ibdmesuc_mayo_score")
sescd_cd  <- parse_characteristics("ibdsescd_totalsescd")

stopifnot(length(titles) == length(gsm),
          length(titles) == length(disease))

# title format: "MSCCR_reGRID_N_Biopsy_M, disease, location inflam tissue"
biopsy_id <- str_extract(titles, "MSCCR_reGRID_\\d+_Biopsy_\\d+")

pdata <- data.frame(
  gsm_id       = gsm,
  biopsy_id    = biopsy_id,
  disease      = disease,
  location     = location,
  inflammation = inflam,
  disease_type = disease2,
  age          = suppressWarnings(as.numeric(age)),
  gender       = gender,
  endo_sev     = endo_sev,
  nancy        = suppressWarnings(as.numeric(nancy)),
  mayo_uc      = suppressWarnings(as.numeric(mayo_uc)),
  sescd_cd     = suppressWarnings(as.numeric(sescd_cd)),
  stringsAsFactors = FALSE
)

cat(sprintf("   샘플 수: %d\n", nrow(pdata)))
cat("   disease 분포:\n"); print(table(pdata$disease))
cat("   location 분포:\n"); print(table(pdata$location))
cat("   inflammation 분포:\n"); print(table(pdata$inflammation))
cat("   disease × inflammation 분포:\n"); print(table(pdata$disease, pdata$inflammation))

write.csv(pdata, "results/00_pdata_clean.csv", row.names = FALSE)
saveRDS(pdata, "rds/pdata.rds")

# ==============================================================================
# 2. Count matrix 로드 (data.table::fread)
# ==============================================================================
cat(">>> [2/8] Loading counts matrix (101MB gz, expect 1-2 min)...\n")
counts_file <- "data/GSE193677_counts.txt.gz"
# 파일 구조: header에 sample ID 2490개, 데이터 행은 "geneID val1 val2 ..." (총 2491토큰)
# fread 자동 감지가 안 정확하므로 header와 body 명시적으로 분리

# 1) Header 읽기 — sample IDs
hdr_conn <- gzfile(counts_file, "r")
hdr_line <- readLines(hdr_conn, n = 1)
close(hdr_conn)
sample_ids <- strsplit(hdr_line, " ")[[1]]
sample_ids <- gsub('^"|"$', '', sample_ids)
cat(sprintf("   header 샘플 수: %d\n", length(sample_ids)))
cat(sprintf("   첫 3개: %s\n", paste(head(sample_ids, 3), collapse = " | ")))

# 2) Body — 첫 행 skip
counts_dt <- fread(counts_file, sep = " ", header = FALSE, skip = 1,
                   data.table = FALSE, check.names = FALSE, quote = '"')
cat(sprintf("   body dim: %d x %d\n", nrow(counts_dt), ncol(counts_dt)))
stopifnot(ncol(counts_dt) == length(sample_ids) + 1)

gene_ids <- counts_dt[, 1]
counts_mat <- as.matrix(counts_dt[, -1])
storage.mode(counts_mat) <- "integer"
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sample_ids
cat(sprintf("   첫 gene: %s, 첫 sample col: %s\n",
            head(gene_ids, 1), head(sample_ids, 1)))

cat(sprintf("   counts: %d genes x %d samples\n", nrow(counts_mat), ncol(counts_mat)))

# ==============================================================================
# 3. Sample matching — counts col names vs pdata biopsy_id
# ==============================================================================
cat(">>> [3/8] Sample matching\n")
common <- intersect(colnames(counts_mat), pdata$biopsy_id)
cat(sprintf("   매칭: %d / %d (counts), %d / %d (pdata)\n",
            length(common), ncol(counts_mat),
            length(common), nrow(pdata)))

# 일부 안 맞는 경우 inspect
if (length(common) < 0.9 * ncol(counts_mat)) {
  cat("   counts에만 있는 ID 예시:\n")
  print(head(setdiff(colnames(counts_mat), pdata$biopsy_id)))
  cat("   pdata에만 있는 ID 예시:\n")
  print(head(setdiff(pdata$biopsy_id, colnames(counts_mat))))
}

counts_mat <- counts_mat[, common]
pdata <- pdata[match(common, pdata$biopsy_id), ]
stopifnot(identical(colnames(counts_mat), pdata$biopsy_id))

# ==============================================================================
# 4. QC + 필터링
# ==============================================================================
cat(">>> [4/8] QC: low-count gene & low-library 제거\n")
keep_g <- rowSums(counts_mat >= 10) >= 0.1 * ncol(counts_mat)
counts_mat <- counts_mat[keep_g, ]
cat(sprintf("   gene filter 후: %d genes\n", nrow(counts_mat)))

lib_size <- colSums(counts_mat)
p_lib <- ggplot(data.frame(lib = lib_size / 1e6), aes(lib)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +
  labs(x = "Library size (M reads)", title = "GSE193677 library sizes") +
  theme_minimal()
ggsave("figures/01_library_size.png", p_lib, width = 7, height = 4, dpi = 120)

keep_s <- lib_size >= 0.2 * median(lib_size)
counts_mat <- counts_mat[, keep_s]; pdata <- pdata[keep_s, ]
cat(sprintf("   sample QC 후: %d samples\n", ncol(counts_mat)))

# ==============================================================================
# 5. DESeq2 vst (single GSE이므로 ComBat-seq skip)
# ==============================================================================
expr_cache <- "rds/expr_vst.rds"
pdata_cache <- "rds/pdata_final.rds"
if (file.exists(expr_cache) && file.exists(pdata_cache)) {
  cat(">>> [5/8] Loading cached vst result (skip)\n")
  expr_mat <- readRDS(expr_cache)
  pdata <- readRDS(pdata_cache)
} else {
  cat(">>> [5/8] DESeq2 normalization + vst (수십 분 소요)\n")
  pdata$disease <- factor(pdata$disease, levels = c("Control", "UC", "CD"))
  pdata$inflammation <- factor(pdata$inflammation)

  dds <- DESeqDataSetFromMatrix(counts_mat, colData = pdata, design = ~ disease)
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = FALSE)
  expr_mat <- assay(vsd)

  # Ensembl → Symbol
  cat("   Ensembl → SYMBOL 변환\n")
  ens_clean <- sub("\\..*$", "", rownames(expr_mat))
  sym_map <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
           column = "SYMBOL", multiVals = "first"))
  new_rn <- ifelse(is.na(sym_map), ens_clean, sym_map)
  v <- matrixStats::rowVars(expr_mat)
  ord <- order(v, decreasing = TRUE)
  expr_mat <- expr_mat[ord, ]; new_rn <- new_rn[ord]
  expr_mat <- expr_mat[!duplicated(new_rn), ]
  rownames(expr_mat) <- new_rn[!duplicated(new_rn)]
  cat(sprintf("   최종 expression matrix: %d genes x %d samples\n",
              nrow(expr_mat), ncol(expr_mat)))

  saveRDS(expr_mat, expr_cache)
  saveRDS(pdata, pdata_cache)
}

# ==============================================================================
# 6. SIRT6 anchor 상관관계
# ==============================================================================
cat(">>> [6/8] SIRT6 anchor Spearman correlation\n")
stopifnot("SIRT6" %in% rownames(expr_mat))

target_genes <- c(NLRP3 = "NLRP3", CASP1 = "CASP1", IL1B = "IL1B",
                  ASC = "PYCARD", TAK1 = "MAP3K7", NFKB_p65 = "RELA",
                  FOXC1 = "FOXC1")

# Subgroup 정의
subgroups <- list(
  All_samples    = pdata$biopsy_id,
  Control        = pdata$biopsy_id[pdata$disease == "Control"],
  UC_all         = pdata$biopsy_id[pdata$disease == "UC"],
  CD_all         = pdata$biopsy_id[pdata$disease == "CD"],
  UC_inflamed    = pdata$biopsy_id[pdata$disease == "UC" & pdata$inflammation == "I"],
  UC_noninflamed = pdata$biopsy_id[pdata$disease == "UC" & pdata$inflammation == "NonI"],
  CD_inflamed    = pdata$biopsy_id[pdata$disease == "CD" & pdata$inflammation == "I"],
  CD_noninflamed = pdata$biopsy_id[pdata$disease == "CD" & pdata$inflammation == "NonI"]
)
cat("   Subgroup sizes:\n")
for (n in names(subgroups)) cat(sprintf("     %s: n=%d\n", n, length(subgroups[[n]])))

compute_cor <- function(samps, label) {
  if (length(samps) < 10) return(NULL)
  e <- expr_mat[, samps]
  s <- e["SIRT6", ]
  out <- lapply(names(target_genes), function(tg) {
    g <- target_genes[tg]
    if (!g %in% rownames(e)) return(NULL)
    ct <- cor.test(s, e[g, ], method = "spearman", exact = FALSE)
    data.frame(subset = label, n = length(samps),
               target_label = tg, target_gene = unname(g),
               rho = unname(ct$estimate), p = ct$p.value)
  })
  bind_rows(out)
}

cor_results <- bind_rows(lapply(names(subgroups), function(n)
  compute_cor(subgroups[[n]], n))) %>%
  group_by(subset) %>%
  mutate(q = p.adjust(p, method = "BH")) %>% ungroup() %>%
  arrange(subset, p)

write.csv(cor_results, "results/01_sirt6_correlation_table.csv", row.names = FALSE)

# Heatmap
cor_wide <- cor_results %>% select(subset, target_label, rho) %>%
  pivot_wider(names_from = target_label, values_from = rho) %>%
  column_to_rownames("subset") %>% as.matrix()
pheatmap(cor_wide,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-0.6, 0.6, length.out = 101),
  main = "SIRT6 vs target genes — Spearman ρ (GSE193677)",
  cluster_rows = FALSE, cluster_cols = TRUE,
  display_numbers = TRUE, fontsize_number = 9,
  filename = "figures/02_sirt6_correlation_heatmap.png",
  width = 8, height = 5)

# UC_inflamed scatter
uc_i <- subgroups$UC_inflamed
if (length(uc_i) >= 30) {
  scatter <- lapply(names(target_genes), function(tg) {
    g <- target_genes[tg]
    df <- data.frame(SIRT6 = expr_mat["SIRT6", uc_i],
                     target = expr_mat[g, uc_i])
    ct <- cor.test(df$SIRT6, df$target, method = "spearman", exact = FALSE)
    ggplot(df, aes(SIRT6, target)) +
      geom_point(alpha = 0.4, size = 0.8, color = "steelblue") +
      geom_smooth(method = "lm", color = "darkred") +
      labs(x = "SIRT6 (vst)", y = sprintf("%s (%s)", tg, g),
           title = sprintf("ρ=%.3f, q=%.2e", ct$estimate, ct$p.value)) +
      theme_minimal(base_size = 9)
  })
  ggsave("figures/03_sirt6_scatter_UC_inflamed.png",
         ggarrange(plotlist = scatter, ncol = 4, nrow = 2),
         width = 14, height = 7, dpi = 120)
}

# ==============================================================================
# 7. WGCNA
# ==============================================================================
cat(">>> [7/8] WGCNA (UC_inflamed + Control subset, 30~90분)\n")

focus_samples <- c(subgroups$UC_inflamed, subgroups$Control)
focus_samples <- intersect(focus_samples, colnames(expr_mat))

# Top variance genes
v <- matrixStats::rowVars(expr_mat[, focus_samples])
var_idx <- order(v, decreasing = TRUE)[1:5000]
sirt6_idx <- which(rownames(expr_mat) == "SIRT6")
var_idx <- unique(c(var_idx, sirt6_idx))
expr_w <- t(expr_mat[var_idx, focus_samples])

cat(sprintf("   WGCNA 입력: %d samples x %d genes\n", nrow(expr_w), ncol(expr_w)))

powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(expr_w, powerVector = powers,
                         networkType = "signed", verbose = 0)
chosen_power <- if (is.na(sft$powerEstimate)) 12 else sft$powerEstimate
cat(sprintf("   power = %d\n", chosen_power))

png("figures/04_wgcna_sft.png", 1200, 500, res = 130)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     type = "b", xlab = "Power", ylab = "Scale-free R^2",
     main = "Scale-free fit"); abline(h = 0.85, col = "red")
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     type = "b", xlab = "Power", ylab = "Mean connectivity",
     main = "Mean connectivity")
dev.off()

net <- blockwiseModules(expr_w, power = chosen_power,
                        networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, maxBlockSize = 6000,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        verbose = 0, saveTOMs = FALSE)
module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_w)
sirt6_module_color <- unname(module_colors["SIRT6"])
sirt6_module_size  <- sum(module_colors == sirt6_module_color)
cat(sprintf("   SIRT6 모듈: %s (size = %d)\n",
            sirt6_module_color, sirt6_module_size))

# Hub gene
sirt6_module_label <- net$colors["SIRT6"]
mod_genes <- names(net$colors)[net$colors == sirt6_module_label]
ME <- net$MEs[, paste0("ME", sirt6_module_label)]
kME <- apply(expr_w[, mod_genes], 2, function(x) cor(x, ME))
hub <- data.frame(gene = names(kME), kME = kME) %>%
  arrange(desc(abs(kME))) %>% head(30)
write.csv(hub, "results/02_sirt6_module_hub_top30.csv", row.names = FALSE)
cat("   Top 30 hub gene:\n"); print(hub)

# Enrichment
cat(">>> [8/8] Pathway enrichment\n")
ent_ids <- mapIds(org.Hs.eg.db, keys = mod_genes, keytype = "SYMBOL",
                  column = "ENTREZID", multiVals = "first") %>% na.omit() %>% unique()
go_bp <- enrichGO(ent_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
kegg  <- enrichKEGG(ent_ids, organism = "hsa",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.1)
write.csv(as.data.frame(go_bp@result),
          "results/03_GO_BP_sirt6_module.csv", row.names = FALSE)
write.csv(as.data.frame(kegg@result),
          "results/04_KEGG_sirt6_module.csv", row.names = FALSE)

if (nrow(go_bp@result) > 0) {
  png("figures/05_GO_BP_top10.png", 1400, 700, res = 130)
  print(barplot(go_bp, showCategory = 10, title = "GO BP — SIRT6 module"))
  dev.off()
}
if (nrow(kegg@result) > 0) {
  png("figures/06_KEGG_top10.png", 1400, 700, res = 130)
  print(barplot(kegg, showCategory = 10, title = "KEGG — SIRT6 module"))
  dev.off()
}

# ==============================================================================
# Go/No-go 자동 판정
# ==============================================================================
THRESH <- list(rho_strong = 0.25, rho_modest = 0.15,
               q_cutoff = 0.05, q_cutoff_loose = 0.10,
               go_n_strong = 3)

uc_i_res <- cor_results %>% filter(subset == "UC_inflamed")
n_strong <- sum(abs(uc_i_res$rho) >= THRESH$rho_strong &
                uc_i_res$q < THRESH$q_cutoff, na.rm = TRUE)
n_modest_any <- cor_results %>%
  filter(abs(rho) >= THRESH$rho_modest, q < THRESH$q_cutoff_loose) %>% nrow()

wgcna_signal <- sirt6_module_color != "grey" &&
                any(grepl("inflam|immune|cytokine|NF.kB|nlrp|interleu",
                          go_bp@result$Description[go_bp@result$p.adjust < 0.05],
                          ignore.case = TRUE))

decision <- if (n_strong >= THRESH$go_n_strong && wgcna_signal) "GO" else
            if (n_modest_any == 0 && !wgcna_signal) "NO-GO" else
            "CAUTION"

cat("\n======================================\n")
cat(sprintf("UC_inflamed strong corr (|ρ|≥0.25, q<0.05): %d / 7\n", n_strong))
cat(sprintf("Any subgroup modest corr: %d\n", n_modest_any))
cat(sprintf("WGCNA SIRT6 module: %s, inflam enriched: %s\n",
            sirt6_module_color, wgcna_signal))
cat(sprintf("DECISION: %s\n", decision))
cat("======================================\n")

saveRDS(list(
  geo_id = "GSE193677", date = Sys.Date(),
  n_samples_qc = ncol(expr_mat),
  subgroup_n = sapply(subgroups, length),
  cor_results = cor_results,
  hub_top30 = hub,
  sirt6_module = list(color = sirt6_module_color, size = sirt6_module_size),
  go_bp_top = head(go_bp@result, 10),
  kegg_top = head(kegg@result, 10),
  decision = decision,
  decision_inputs = list(n_strong_UC_inflamed = n_strong,
                         n_modest_any = n_modest_any,
                         wgcna_signal = wgcna_signal)
), "rds/pilot_summary.rds")

# Auto-render report
if (file.exists("02_pilot_report.Rmd")) {
  rmarkdown::render("02_pilot_report.Rmd",
                    output_file = "SIRT6_pilot_report.html",
                    output_dir = "results", quiet = TRUE)
}
cat("=== DONE ===\n")
