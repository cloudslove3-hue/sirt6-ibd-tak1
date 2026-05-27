# ==============================================================================
# Master runner — scRNA 분석 전체 파이프라인 순차 실행
# 권장 실행 순서:
#   1. 00_install_scrna_packages.R (한 번)
#   2. 01_smillie_seurat_analysis.R (~30 min, 11K cells × Seurat 5 + SingleR)
#   3. 02_music_deconvolution.R     (~20 min, bulk 2,483 sample × MuSiC)
#   4. 04_cibersortx_signature_export.R (~5 min, CIBERSORTx 업로드용 매트릭스)
#   5. 03_hpa_ihc_download.R        (~5 min, HPA 이미지 다운로드)
#
# 캐싱: Seurat 객체 + MuSiC deconv 모두 rds/에 저장. 재실행 시 빠름.
# ==============================================================================

Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

scripts <- c(
  "01_smillie_seurat_analysis.R",
  "02_music_deconvolution.R",
  "04_cibersortx_signature_export.R",
  "05_rela_target_decoupling.R",
  "03_hpa_ihc_download.R"
)

for (s in scripts) {
  cat(sprintf("\n\n##### Running: %s #####\n\n", s))
  tryCatch(source(s, echo = FALSE),
           error = function(e) {
             cat(sprintf("ERROR in %s: %s\n", s, conditionMessage(e)))
           })
}

cat("\n=== ALL SCRIPTS COMPLETE ===\n")
