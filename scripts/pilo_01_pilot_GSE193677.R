# ==============================================================================
# SIRT6 1차 시범 분석 — GSE193677 (Mount Sinai IBD cohort, n≈2,490)
# 목적: SIRT6 vs {NLRP3, CASP1, IL1B, PYCARD(ASC), MAP3K7(TAK1), RELA, FOXC1}
#       Spearman 상관 + WGCNA 모듈 동정 → go/no-go 판정
# 환경: 로컬 R 4.3+ / RStudio / Windows / 32GB RAM
# 작성: 2026-05-21
# 실행 시간 예상: 다운로드 30분~1h + DESeq2 vst 30~60분 + WGCNA 1~3h
# ==============================================================================

# ---- 0. 작업 디렉토리 & 출력 폴더 ----
# RStudio에서: Session > Set Working Directory > To Source File Location
# 또는 setwd() 수동 지정
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/pilot_analysis")
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("rds", showWarnings = FALSE)

# ---- 1. 패키지 설치 & 로드 ----
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("GEOquery", "DESeq2", "limma", "edgeR", "sva",
               "WGCNA", "clusterProfiler", "org.Hs.eg.db",
               "ReactomePA", "AnnotationDbi", "biomaRt")
cran_pkgs <- c("tidyverse", "data.table", "ggplot2", "ggpubr", "corrplot",
               "pheatmap", "RColorBrewer", "ggrepel", "rmarkdown", "knitr")

for (p in bioc_pkgs) if (!require(p, character.only = TRUE, quietly = TRUE))
  BiocManager::install(p, update = FALSE, ask = FALSE)
for (p in cran_pkgs) if (!require(p, character.only = TRUE, quietly = TRUE))
  install.packages(p, dependencies = TRUE)

suppressPackageStartupMessages({
  library(GEOquery); library(DESeq2); library(limma); library(edgeR)
  library(sva); library(WGCNA); library(clusterProfiler); library(org.Hs.eg.db)
  library(ReactomePA); library(tidyverse); library(data.table)
  library(ggpubr); library(corrplot); library(pheatmap); library(ggrepel)
})

# WGCNA 권장 설정
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))

# ==============================================================================
# 2. 데이터 다운로드 — GSE193677
# ==============================================================================
# GSE193677은 raw/adjusted counts를 supplementary file로 제공.
# GEOquery의 getGEO()는 series matrix(메타데이터+normalized expression)만 받음.
# Raw count matrix는 별도 supplementary file download 필요.

GEO_ID <- "GSE193677"
geo_dir <- "data"

# 2-1. 메타데이터 (series matrix)
message(">>> [1/8] Series matrix 다운로드 — pData 추출용")
if (!file.exists(file.path(geo_dir, "geo_series.rds"))) {
  gse <- getGEO(GEO_ID, destdir = geo_dir, GSEMatrix = TRUE,
                getGPL = FALSE)
  saveRDS(gse, file.path(geo_dir, "geo_series.rds"))
} else {
  gse <- readRDS(file.path(geo_dir, "geo_series.rds"))
}

# 일반적으로 list of ExpressionSet — 첫 element 사용
eset <- gse[[1]]
pdata_raw <- pData(eset)

message(sprintf("   샘플 수: %d, 메타 컬럼 수: %d",
                nrow(pdata_raw), ncol(pdata_raw)))
write.csv(pdata_raw, "data/pdata_raw.csv", row.names = FALSE)

# 2-2. 메타데이터 정제 — 사용자가 GEO에서 직접 컬럼명 확인 필요할 수 있음
# 아래 컬럼 매핑은 일반적인 GSE193677 구조를 가정. 실제와 다르면 조정.
message(">>> [2/8] 메타데이터 정제")

# GEO characteristics 컬럼은 보통 "characteristics_ch1.x" 또는 자유 텍스트 컬럼
# 흔히 사용되는 매핑 시도 — 실제 컬럼명은 names(pdata_raw)로 확인 후 수정
candidate_disease_cols <- grep("disease|diagnosis|status",
                               names(pdata_raw), ignore.case = TRUE, value = TRUE)
candidate_location_cols <- grep("location|tissue|site|region|biopsy",
                                names(pdata_raw), ignore.case = TRUE, value = TRUE)
candidate_inflam_cols <- grep("inflam|active|severity",
                              names(pdata_raw), ignore.case = TRUE, value = TRUE)

message("발견된 후보 컬럼 (실제 데이터 보고 final 매핑 결정):")
message("   disease: ", paste(candidate_disease_cols, collapse = ", "))
message("   location: ", paste(candidate_location_cols, collapse = ", "))
message("   inflammation: ", paste(candidate_inflam_cols, collapse = ", "))

# *** 사용자 직접 확인 단계 ***
# 위 메시지 보고 실제 컬럼명을 아래 변수에 할당 후 진행
COL_DISEASE   <- candidate_disease_cols[1]   # 예: "disease:ch1" 또는 "characteristics_ch1.1"
COL_LOCATION  <- candidate_location_cols[1]
COL_INFLAM    <- candidate_inflam_cols[1]
COL_SAMPLE_ID <- "geo_accession"  # GSM ID

stopifnot(!is.na(COL_DISEASE), !is.na(COL_LOCATION), !is.na(COL_INFLAM))

pdata <- pdata_raw %>%
  as.data.frame() %>%
  rownames_to_column("gsm_id") %>%
  mutate(
    disease      = factor(.data[[COL_DISEASE]]),
    location     = factor(.data[[COL_LOCATION]]),
    inflammation = factor(.data[[COL_INFLAM]])
  ) %>%
  select(gsm_id, disease, location, inflammation, everything())

message(sprintf("   정제 후 샘플 수: %d", nrow(pdata)))
message("   disease 분포:")
print(table(pdata$disease))
message("   location 분포:")
print(table(pdata$location))
message("   inflammation 분포:")
print(table(pdata$inflammation))

saveRDS(pdata, "rds/pdata_cleaned.rds")

# 2-3. Raw count matrix 다운로드 (supplementary)
message(">>> [3/8] Raw count matrix 다운로드 (수 GB일 수 있음)")
sup_dir <- file.path(geo_dir, GEO_ID)
if (!dir.exists(sup_dir)) {
  getGEOSuppFiles(GEO_ID, baseDir = geo_dir, makeDirectory = TRUE)
}
list.files(sup_dir)

# 일반적으로 .txt.gz 또는 .csv.gz 형태. 파일명에 "raw" 또는 "counts" 포함.
counts_file <- list.files(sup_dir, pattern = "(raw|count)", full.names = TRUE,
                          ignore.case = TRUE)[1]
message("   사용할 count file: ", basename(counts_file))

# fread + automagic 압축 처리
counts_mat <- fread(counts_file, data.table = FALSE)
# 첫 컬럼은 보통 gene id (Ensembl 또는 symbol)
rownames(counts_mat) <- counts_mat[, 1]
counts_mat <- as.matrix(counts_mat[, -1])
storage.mode(counts_mat) <- "integer"  # DESeq2는 integer 요구
message(sprintf("   counts 차원: %d genes x %d samples",
                nrow(counts_mat), ncol(counts_mat)))

# Sample ID 매칭
common_samples <- intersect(colnames(counts_mat), pdata$gsm_id)
if (length(common_samples) < 0.8 * ncol(counts_mat)) {
  warning("샘플 ID 매칭률이 낮음. GSE193677 supplementary는 column명이 GSM이 아닌
           다른 ID(예: subject_ID)일 수 있음. pdata의 다른 컬럼과 매칭 시도 필요.")
  # fallback: GSE193677은 종종 'title' 컬럼에 sample명이 있음
  matched <- match(colnames(counts_mat), pdata$title)
  if (sum(!is.na(matched)) > 0.8 * ncol(counts_mat)) {
    pdata$gsm_id <- pdata$title[matched]
    common_samples <- intersect(colnames(counts_mat), pdata$gsm_id)
  }
}
message(sprintf("   매칭된 샘플: %d", length(common_samples)))

counts_mat <- counts_mat[, common_samples]
pdata <- pdata[match(common_samples, pdata$gsm_id), ]
stopifnot(identical(colnames(counts_mat), pdata$gsm_id))

saveRDS(counts_mat, "rds/counts_raw.rds")

# ==============================================================================
# 3. QC 및 필터링
# ==============================================================================
message(">>> [4/8] QC: low-count gene 제거, library size 확인")

# 최소 expression filter
keep <- rowSums(counts_mat >= 10) >= 0.1 * ncol(counts_mat)
counts_mat <- counts_mat[keep, ]
message(sprintf("   filtering 후: %d genes", nrow(counts_mat)))

# Library size 분포
lib_size <- colSums(counts_mat)
p_lib <- ggplot(data.frame(lib_size = lib_size), aes(x = lib_size / 1e6)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(x = "Library size (M reads)", y = "Sample count",
       title = "Library size distribution") +
  theme_minimal()
ggsave("figures/01_library_size.png", p_lib, width = 7, height = 4, dpi = 150)

# 극단적으로 낮은 library 제거 (median의 20% 미만)
keep_sample <- lib_size >= 0.2 * median(lib_size)
counts_mat <- counts_mat[, keep_sample]
pdata <- pdata[keep_sample, ]
message(sprintf("   QC pass 샘플 수: %d", ncol(counts_mat)))

# ==============================================================================
# 4. DESeq2 정규화 + vst
# ==============================================================================
message(">>> [5/8] DESeq2 정규화 + variance stabilizing transform")

# ComBat-seq 자리 — 단일 GSE에서는 적용하지 않음.
# 다음에 GSE99816 등을 합칠 때 아래 함수 호출:
#   counts_adj <- sva::ComBat_seq(counts = combined_counts,
#                                  batch = combined_meta$dataset,
#                                  group = combined_meta$disease)
# 이번 시범 분석에서는 skip.

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = pdata,
  design    = ~ disease  # 다른 covariate은 sub-analysis에서 추가
)
dds <- estimateSizeFactors(dds)

# vst — variance stabilizing, WGCNA & 상관분석에 권장
vsd <- vst(dds, blind = FALSE)
expr_mat <- assay(vsd)  # gene x sample, log2-scale 유사
saveRDS(vsd, "rds/vsd.rds")
saveRDS(expr_mat, "rds/expr_vst.rds")

# Gene ID 변환 — Ensembl이면 symbol로
if (any(grepl("^ENSG", rownames(expr_mat)))) {
  message("   Ensembl ID 감지 — symbol로 변환")
  ens_clean <- sub("\\..*$", "", rownames(expr_mat))
  symbol_map <- mapIds(org.Hs.eg.db, keys = ens_clean,
                       keytype = "ENSEMBL", column = "SYMBOL",
                       multiVals = "first")
  rownames(expr_mat) <- ifelse(is.na(symbol_map), ens_clean, symbol_map)
  # 중복 제거 (가장 발현 높은 것 유지)
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)) & !is.na(rownames(expr_mat)), ]
}

message(sprintf("   최종 expression matrix: %d genes x %d samples",
                nrow(expr_mat), ncol(expr_mat)))

# ==============================================================================
# 5. SIRT6 anchor 상관관계 분석
# ==============================================================================
message(">>> [6/8] SIRT6 anchor correlation")

if (!"SIRT6" %in% rownames(expr_mat)) {
  stop("SIRT6가 expression matrix에 없습니다. ID 매핑 재확인 필요.")
}

target_genes <- c(
  "NLRP3"  = "NLRP3",
  "CASP1"  = "CASP1",
  "IL1B"   = "IL1B",
  "ASC"    = "PYCARD",
  "TAK1"   = "MAP3K7",
  "NFKB_p65" = "RELA",
  "FOXC1"  = "FOXC1"
)

# 함수 — subset별 Spearman 상관
compute_sirt6_cor <- function(expr, samples, label) {
  e <- expr[, samples, drop = FALSE]
  s <- e["SIRT6", ]
  out <- lapply(names(target_genes), function(tg) {
    g <- target_genes[[tg]]
    if (!g %in% rownames(e)) return(NULL)
    ct <- cor.test(s, e[g, ], method = "spearman", exact = FALSE)
    data.frame(
      subset = label, n = length(samples),
      target_label = tg, target_gene = g,
      rho = unname(ct$estimate), p = ct$p.value
    )
  })
  bind_rows(out)
}

# 5-1. 전체 + 주요 subgroup
subgroups <- list(
  "All_samples"    = pdata$gsm_id,
  "Control"        = pdata$gsm_id[grepl("control|ctrl|healthy|non.?IBD",
                                         pdata$disease, ignore.case = TRUE)],
  "UC_all"         = pdata$gsm_id[grepl("UC|ulcerative", pdata$disease, ignore.case = TRUE)],
  "CD_all"         = pdata$gsm_id[grepl("CD|crohn", pdata$disease, ignore.case = TRUE)],
  "UC_inflamed"    = pdata$gsm_id[grepl("UC|ulcerative", pdata$disease, ignore.case = TRUE) &
                                   grepl("^I$|inflamed", pdata$inflammation, ignore.case = TRUE) &
                                   !grepl("non", pdata$inflammation, ignore.case = TRUE)],
  "UC_noninflamed" = pdata$gsm_id[grepl("UC|ulcerative", pdata$disease, ignore.case = TRUE) &
                                   grepl("non.?inflamed|NonI", pdata$inflammation, ignore.case = TRUE)],
  "CD_inflamed"    = pdata$gsm_id[grepl("CD|crohn", pdata$disease, ignore.case = TRUE) &
                                   grepl("^I$|inflamed", pdata$inflammation, ignore.case = TRUE) &
                                   !grepl("non", pdata$inflammation, ignore.case = TRUE)],
  "CD_noninflamed" = pdata$gsm_id[grepl("CD|crohn", pdata$disease, ignore.case = TRUE) &
                                   grepl("non.?inflamed|NonI", pdata$inflammation, ignore.case = TRUE)]
)
# 각 subgroup 샘플 수 출력
for (n in names(subgroups))
  message(sprintf("   %s: n = %d", n, length(subgroups[[n]])))

cor_results <- bind_rows(lapply(names(subgroups), function(label) {
  s <- subgroups[[label]]
  if (length(s) < 10) return(NULL)
  compute_sirt6_cor(expr_mat, s, label)
})) %>%
  group_by(subset) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(subset, p)

write.csv(cor_results, "results/01_sirt6_correlation_table.csv", row.names = FALSE)
saveRDS(cor_results, "rds/cor_results.rds")

# 5-2. 시각화 — heatmap (subset × target gene)
cor_wide <- cor_results %>%
  select(subset, target_label, rho) %>%
  pivot_wider(names_from = target_label, values_from = rho) %>%
  column_to_rownames("subset") %>%
  as.matrix()

p_cor_heat <- pheatmap(cor_wide,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-0.6, 0.6, length.out = 101),
  main = "SIRT6 vs target genes — Spearman ρ",
  cluster_rows = FALSE, cluster_cols = TRUE,
  display_numbers = TRUE, fontsize_number = 8,
  filename = "figures/02_sirt6_correlation_heatmap.png",
  width = 8, height = 5
)

# 5-3. UC_inflamed subset scatter plot 7개
uc_infl_samples <- subgroups$UC_inflamed
if (length(uc_infl_samples) >= 30) {
  s <- expr_mat["SIRT6", uc_infl_samples]
  scatter_plots <- lapply(names(target_genes), function(tg) {
    g <- target_genes[[tg]]
    if (!g %in% rownames(expr_mat)) return(NULL)
    df <- data.frame(SIRT6 = s, target = expr_mat[g, uc_infl_samples])
    ct <- cor.test(df$SIRT6, df$target, method = "spearman", exact = FALSE)
    ggplot(df, aes(SIRT6, target)) +
      geom_point(alpha = 0.5, color = "steelblue", size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "darkred") +
      labs(x = "SIRT6 (vst)", y = paste0(tg, " (", g, ", vst)"),
           title = sprintf("ρ = %.3f, p = %.2e", ct$estimate, ct$p.value)) +
      theme_minimal(base_size = 9)
  })
  scatter_plots <- scatter_plots[!sapply(scatter_plots, is.null)]
  g_combined <- ggarrange(plotlist = scatter_plots, ncol = 4, nrow = 2)
  ggsave("figures/03_sirt6_scatter_UC_inflamed.png", g_combined,
         width = 14, height = 7, dpi = 150)
}

# ==============================================================================
# 6. WGCNA — SIRT6 module 동정
# ==============================================================================
message(">>> [7/8] WGCNA (수십분~수시간 소요)")

# WGCNA는 sample x gene 매트릭스 요구
# 메모리 효율을 위해 variance 상위 5,000개 gene만 사용 (WGCNA 권장)
var_genes <- order(apply(expr_mat, 1, var), decreasing = TRUE)[1:5000]
# SIRT6는 반드시 포함
sirt6_idx <- which(rownames(expr_mat) == "SIRT6")
var_genes <- unique(c(var_genes, sirt6_idx))
expr_wgcna <- t(expr_mat[var_genes, ])  # sample x gene

# 6-1. 1차 분석 코호트 선택 — 신호가 가장 명확할 UC_inflamed + Control
focus_samples <- c(subgroups$UC_inflamed, subgroups$Control)
focus_samples <- intersect(focus_samples, rownames(expr_wgcna))
expr_wgcna_focus <- expr_wgcna[focus_samples, ]
message(sprintf("   WGCNA 입력: %d samples x %d genes",
                nrow(expr_wgcna_focus), ncol(expr_wgcna_focus)))

# 6-2. soft-thresholding power 선택
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(expr_wgcna_focus, powerVector = powers,
                         networkType = "signed", verbose = 0)
# scale-free R^2 ≥ 0.85 만족하는 최소 power
chosen_power <- sft$powerEstimate
if (is.na(chosen_power)) chosen_power <- 12
message(sprintf("   선택된 soft-thresholding power: %d", chosen_power))

png("figures/04_wgcna_sft.png", width = 1000, height = 500, res = 120)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power,
     -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold (power)", ylab = "Scale Free R^2",
     type = "b", main = "Scale-free fit")
abline(h = 0.85, col = "red")
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)", ylab = "Mean connectivity",
     type = "b", main = "Mean connectivity")
dev.off()

# 6-3. blockwise module detection (메모리 절약)
net <- blockwiseModules(
  expr_wgcna_focus,
  power = chosen_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  maxBlockSize = 5000,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "rds/TOM",
  verbose = 3
)
saveRDS(net, "rds/wgcna_net.rds")

module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(expr_wgcna_focus)
sirt6_module_color <- module_colors["SIRT6"]
message(sprintf("   SIRT6가 속한 module: %s (size = %d)",
                sirt6_module_color, sum(module_colors == sirt6_module_color)))

# 6-4. SIRT6 module의 hub gene 30개 (intramodular connectivity)
ME <- net$MEs
sirt6_module_label <- net$colors["SIRT6"]
sirt6_module_genes <- names(net$colors)[net$colors == sirt6_module_label]
sirt6_module_expr <- expr_wgcna_focus[, sirt6_module_genes]

# kME = correlation with module eigengene
ME_sirt6 <- ME[, paste0("ME", sirt6_module_label)]
kME <- apply(sirt6_module_expr, 2, function(x) cor(x, ME_sirt6, method = "pearson"))
hub_df <- data.frame(gene = names(kME), kME = kME) %>%
  arrange(desc(abs(kME)))
top30_hub <- head(hub_df, 30)
write.csv(top30_hub, "results/02_sirt6_module_hub_top30.csv", row.names = FALSE)

message("   SIRT6 module hub gene (top 30):")
print(top30_hub)

# 6-5. GO/KEGG/Reactome enrichment
message(">>> [8/8] Pathway enrichment (SIRT6 module)")

ent_ids <- mapIds(org.Hs.eg.db, keys = sirt6_module_genes,
                  keytype = "SYMBOL", column = "ENTREZID",
                  multiVals = "first") %>% na.omit()

go_bp <- enrichGO(gene = ent_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
kegg <- enrichKEGG(gene = ent_ids, organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.1)
reactome <- enrichPathway(gene = ent_ids, pvalueCutoff = 0.05, readable = TRUE)

write.csv(as.data.frame(go_bp@result),    "results/03_GO_BP_sirt6_module.csv", row.names = FALSE)
write.csv(as.data.frame(kegg@result),     "results/04_KEGG_sirt6_module.csv",  row.names = FALSE)
write.csv(as.data.frame(reactome@result), "results/05_Reactome_sirt6_module.csv", row.names = FALSE)

# 시각화 — top 10
png("figures/05_GO_BP_top10.png", width = 1200, height = 600, res = 120)
print(barplot(go_bp, showCategory = 10, title = "GO BP — SIRT6 module"))
dev.off()

png("figures/06_KEGG_top10.png", width = 1200, height = 600, res = 120)
print(barplot(kegg, showCategory = 10, title = "KEGG — SIRT6 module"))
dev.off()

# ==============================================================================
# 7. Go/No-go 자동 판정
# ==============================================================================
message(">>> Final: Go/No-go 자동 판정")

# 사전 정의된 기준 (Part C 문서와 동일)
GO_NOGO_THRESH <- list(
  rho_strong      = 0.25,   # |ρ| ≥ 0.25 = strong
  rho_modest      = 0.15,   # |ρ| ≥ 0.15 = modest
  q_cutoff        = 0.05,
  q_cutoff_loose  = 0.10,
  go_n_strong     = 3,      # UC_inflamed에서 3개 이상 strong → GO
  nogo_n_modest   = 1       # 어느 subgroup에서도 modest 없으면 → NOGO
)

uc_infl_res <- cor_results %>% filter(subset == "UC_inflamed")
n_strong <- sum(abs(uc_infl_res$rho) >= GO_NOGO_THRESH$rho_strong &
                uc_infl_res$q < GO_NOGO_THRESH$q_cutoff, na.rm = TRUE)
n_modest_any <- cor_results %>%
  filter(abs(rho) >= GO_NOGO_THRESH$rho_modest,
         q < GO_NOGO_THRESH$q_cutoff_loose) %>% nrow()

wgcna_signal <- sirt6_module_color != "grey" &&
                nrow(go_bp@result %>% filter(p.adjust < 0.05 &
                  grepl("inflam|immune|cytokine|NF.kB|nlrp",
                        Description, ignore.case = TRUE))) > 0

decision <- if (n_strong >= GO_NOGO_THRESH$go_n_strong && wgcna_signal) {
  "GO"
} else if (n_modest_any < GO_NOGO_THRESH$nogo_n_modest && !wgcna_signal) {
  "NO-GO"
} else {
  "CAUTION (재평가)"
}

message("======================================")
message(sprintf("UC_inflamed strong correlation (|ρ|≥0.25, q<0.05): %d / 7", n_strong))
message(sprintf("Any subgroup modest correlation (|ρ|≥0.15, q<0.10): %d", n_modest_any))
message(sprintf("WGCNA SIRT6 module: %s, inflammation pathway enriched: %s",
                sirt6_module_color, wgcna_signal))
message(sprintf("DECISION: %s", decision))
message("======================================")

# 결과 요약 저장
summary_obj <- list(
  geo_id = GEO_ID,
  date = Sys.Date(),
  n_samples_qc = ncol(expr_mat),
  subgroup_n = sapply(subgroups, length),
  cor_results = cor_results,
  hub_top30 = top30_hub,
  sirt6_module = list(color = sirt6_module_color,
                      size = sum(module_colors == sirt6_module_color),
                      eigengene_var = var(ME_sirt6)),
  go_bp_top = head(go_bp@result, 10),
  kegg_top = head(kegg@result, 10),
  reactome_top = head(reactome@result, 10),
  decision = decision,
  decision_inputs = list(
    n_strong_UC_inflamed = n_strong,
    n_modest_any = n_modest_any,
    wgcna_signal = wgcna_signal
  )
)
saveRDS(summary_obj, "rds/pilot_summary.rds")

# ==============================================================================
# 8. Auto-render report
# ==============================================================================
message(">>> Rendering 1-page report...")
if (file.exists("02_pilot_report.Rmd")) {
  rmarkdown::render("02_pilot_report.Rmd",
                    output_file = "results/SIRT6_pilot_report.html",
                    quiet = FALSE)
  message("   완성: results/SIRT6_pilot_report.html")
}

message("=== PIPELINE 완료 ===")
message("산출물 위치:")
message("  - results/01_sirt6_correlation_table.csv")
message("  - results/02_sirt6_module_hub_top30.csv")
message("  - results/03_GO_BP_sirt6_module.csv")
message("  - results/04_KEGG_sirt6_module.csv")
message("  - results/05_Reactome_sirt6_module.csv")
message("  - figures/01~06.png")
message("  - results/SIRT6_pilot_report.html")
