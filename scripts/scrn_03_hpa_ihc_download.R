# ==============================================================================
# Human Protein Atlas (HPA) SIRT6 IHC 다운로드 — v2 (tissue-specific pages)
# ==============================================================================
Sys.setenv(TMP = "C:/Temp/R", TEMP = "C:/Temp/R", TMPDIR = "C:/Temp/R")
.libPaths(c("C:/Rlibs", .libPaths()))
setwd("C:/Users/방창석/Dropbox/방창석 2026/33. SIRT6/scrna_analysis")

dir.create("hpa", showWarnings = FALSE)
dir.create("hpa/images", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(httr); library(xml2); library(rvest); library(dplyr); library(stringr)
})

SIRT6_ENSG <- "ENSG00000077463"
tissues <- c("colon", "rectum", "small-intestine", "duodenum",
             "pancreas", "liver")
# (선택) stomach, appendix도 추가 가능

download_log <- data.frame()
all_urls <- list()

for (tissue in tissues) {
  url <- sprintf("https://www.proteinatlas.org/%s-SIRT6/tissue/%s",
                 SIRT6_ENSG, tissue)
  cat(sprintf(">>> %s\n   page: %s\n", tissue, url))

  page <- tryCatch(read_html(url), error = function(e) NULL)
  if (is.null(page)) {
    cat("   FAILED to fetch page\n"); next
  }

  # img tag src 추출
  img_src <- page %>% html_nodes("img") %>% html_attr("src")
  img_src <- img_src[grepl("images\\.proteinatlas\\.org", img_src)]

  # data-* attribute에 있을 수도 있음 (lazy-load)
  data_src <- page %>% html_nodes("[data-src]") %>% html_attr("data-src")
  data_src <- data_src[grepl("images\\.proteinatlas\\.org", data_src)]

  # 모두 합치고 중복 제거
  candidates <- unique(c(img_src, data_src))
  # protocol-relative → absolute
  candidates <- ifelse(startsWith(candidates, "//"),
                       paste0("https:", candidates), candidates)
  # full resolution 우선: '_medium.jpg' 제거
  full_res <- sub("_medium\\.jpg$", ".jpg", candidates)
  candidates <- unique(c(candidates, full_res))

  cat(sprintf("   discovered %d image URLs\n", length(candidates)))
  all_urls[[tissue]] <- candidates

  # 다운로드
  for (u in candidates) {
    # 파일명에 tissue prefix
    fname <- paste0(tissue, "_", basename(u))
    dst <- file.path("hpa/images", fname)
    if (file.exists(dst)) {
      download_log <- rbind(download_log, data.frame(
        tissue = tissue, url = u, file = dst, status = "exists"))
      next
    }
    tryCatch({
      download.file(u, dst, mode = "wb", quiet = TRUE)
      download_log <- rbind(download_log, data.frame(
        tissue = tissue, url = u, file = dst, status = "ok"))
    }, error = function(e) {
      download_log <<- rbind(download_log, data.frame(
        tissue = tissue, url = u, file = dst, status = paste("err:", e$message)))
    })
  }
}

write.csv(download_log, "hpa/download_log.csv", row.names = FALSE)
saveRDS(all_urls, "hpa/url_list.rds")

# 요약
status_summary <- table(download_log$status)
cat("\n=== Download summary ===\n")
print(status_summary)
cat("\nImages saved to: hpa/images/\n")
files <- list.files("hpa/images", pattern = "\\.jpg$")
cat(sprintf("Total %d jpg files\n", length(files)))
if (length(files) > 0) {
  by_tissue <- table(sub("_.*", "", files))
  cat("By tissue:\n"); print(by_tissue)
}

# ==============================================================================
# Single-cell RNA expression from HPA (orthogonal data)
# ==============================================================================
cat("\n>>> HPA single-cell RNA (orthogonal to Smillie scRNA)\n")
hpa_sc_data_url <- "https://www.proteinatlas.org/api/search_download.php?search=gene:SIRT6&format=tsv&columns=g,gs,scal,sccel,sct"
# 사용 가능 컬럼: scal=scRNA average level, sccel=scRNA cell type enriched, sct=scRNA tissue distribution
tryCatch({
  download.file(hpa_sc_data_url, "hpa/SIRT6_singlecell_metadata.tsv",
                mode = "wb", quiet = TRUE)
  cat("   HPA scRNA metadata saved\n")
}, error = function(e) cat("   HPA API fetch failed:", e$message, "\n"))

cat("\n=== HPA IHC v2 DONE ===\n")
cat("Files location:\n")
cat("  hpa/images/*.jpg — tissue-specific IHC staining images\n")
cat("  hpa/download_log.csv — 다운로드 결과\n")
cat("  hpa/url_list.rds — tissue별 image URL 목록\n")
cat("  hpa/SIRT6_singlecell_metadata.tsv — HPA scRNA cell type 정보\n")
