# ==============================================================================
# 脚本名称: 09_build_submission_package.R
# 核心任务: 将所有投稿所需手稿、6 张最终主图、附表 Excel、标准 GMT 基因集与
#           可复现代码统一复制到新建的独立投稿文件夹 06_Journal_Submission_Package/
# 安全承诺: 全程仅使用只读读取与 file.copy()，绝不修改、覆盖或移动任何原有文件！
# ==============================================================================

rm(list = ls())

if (!dir.exists("01_clean_data") && dir.exists("~/HCC_Admixture_Project/01_clean_data")) {
  setwd("~/HCC_Admixture_Project")
}
stopifnot("请先将工作目录切换到 HCC_Admixture_Project 根目录" = dir.exists("01_clean_data"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(tools)
})

# ------------------------------------------------------------------------------
# 1. 创建全新独立的投稿文件夹（若已有同名目录则自动附加时间戳，绝不覆盖旧内容）
# ------------------------------------------------------------------------------
pkg_root <- "06_Journal_Submission_Package"
if (dir.exists(pkg_root)) {
  pkg_root <- paste0(pkg_root, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}

dir_ms      <- file.path(pkg_root, "01_Manuscript")
dir_fig_tif <- file.path(pkg_root, "02_Main_Figures", "TIFF_300dpi")
dir_fig_pdf <- file.path(pkg_root, "02_Main_Figures", "PDF_Vector")
dir_fig_png <- file.path(pkg_root, "02_Main_Figures", "PNG_Preview")
dir_sup     <- file.path(pkg_root, "03_Supplementary_Materials")
dir_sup_csv <- file.path(pkg_root, "03_Supplementary_Materials", "Tables_CSV")
dir_code    <- file.path(pkg_root, "04_Code_and_Reproducibility", "scripts")
dir_repro   <- file.path(pkg_root, "04_Code_and_Reproducibility")

for (d in c(dir_ms, dir_fig_tif, dir_fig_pdf, dir_fig_png, dir_sup, dir_sup_csv, dir_code)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message(">>> [1/6] 已创建全新投稿打包目录: ", pkg_root)

# 记录所有复制/生成日志
manifest_rows <- list()
log_item <- function(category, dest_file, source_info, status) {
  manifest_rows[[length(manifest_rows) + 1]] <<- tibble(
    Category = category,
    Package_Path = dest_file,
    Source = source_info,
    Status = status
  )
}

copy_safe <- function(src, dst, category) {
  if (!is.na(src) && file.exists(src)) {
    ok <- file.copy(from = src, to = dst, overwrite = TRUE, copy.date = TRUE)
    log_item(category, dst, src, ifelse(ok, "COPIED", "COPY_FAILED"))
    return(ok)
  } else {
    log_item(category, dst, ifelse(is.na(src), "NOT_FOUND", src), "MISSING")
    return(FALSE)
  }
}

# ------------------------------------------------------------------------------
# 2. 复制英文终稿 Word (.docx) —— 自动在常见目录检索（只复制，绝不移动）
# ------------------------------------------------------------------------------
message(">>> [2/6] 正在检索并复制最终英文手稿 (.docx)...")
ms_target_name <- "HCC_Admixture_Manuscript_Final_v4_31Refs.docx"
ms_candidates <- c(
  file.path("05_muscript", ms_target_name),
  file.path(".", ms_target_name),
  file.path(path.expand("~/Downloads"), ms_target_name),
  file.path(path.expand("~/Desktop"), ms_target_name)
)
ms_found <- ms_candidates[file.exists(ms_candidates)][1]

if (!is.na(ms_found)) {
  copy_safe(ms_found, file.path(dir_ms, ms_target_name), "Manuscript")
  message("    已找到并复制手稿: ", ms_found)
} else {
  log_item("Manuscript", file.path(dir_ms, ms_target_name),
           "Not found in 05_muscript / root / ~/Downloads / ~/Desktop", "ACTION_NEEDED")
  message("    [提示] 暂未在本地找到 ", ms_target_name,
          "；请下载后手动复制一份放入 ", dir_ms, "/ 即可。")
}

# ------------------------------------------------------------------------------
# 3. 复制 6 张最终主图 (TIFF + PDF + PNG)
#    Figure 1–3: 优先选用 05_muscript/review_20261123/figures 中性措辞版，其次 final_manuscript
#    Figure 4–6: 严格选用 04_figures/final_manuscript_v2 最终终审版
# ------------------------------------------------------------------------------
message(">>> [3/6] 正在复制 6 张最终主图 (TIFF / PDF / PNG)...")

fig_specs <- list(
  list(id = "Figure_1",
       dirs = c("05_muscript/review_20261123/figures", "04_figures/final_manuscript"),
       stem = "Figure1_Admixture_Sensitivity_Landscape"),
  list(id = "Figure_2",
       dirs = c("05_muscript/review_20261123/figures", "04_figures/final_manuscript"),
       stem = "Figure2_Source_and_MultiCohort_Validation"),
  list(id = "Figure_3",
       dirs = c("05_muscript/review_20261123/figures", "04_figures/final_manuscript"),
       stem = "Figure3_scRNA_Mechanism_Decomposition"),
  list(id = "Figure_4",
       dirs = c("04_figures/final_manuscript_v2"),
       stem = "Figure4_Spatial_All7_Validation"),
  list(id = "Figure_5",
       dirs = c("04_figures/final_manuscript_v2"),
       stem = "Figure5_Clinical_Survival_and_NHAS_Diagnostics"),
  list(id = "Figure_6",
       dirs = c("04_figures/final_manuscript_v2"),
       stem = "Figure6_Benchmark_Null_and_Multivar_Survival")
)

ext_dest_map <- list(
  tiff = dir_fig_tif,
  pdf  = dir_fig_pdf,
  png  = dir_fig_png
)

for (fs in fig_specs) {
  picked_dir <- NA_character_
  for (cand_d in fs$dirs) {
    if (file.exists(file.path(cand_d, paste0(fs$stem, ".pdf"))) ||
        file.exists(file.path(cand_d, paste0(fs$stem, ".png")))) {
      picked_dir <- cand_d
      break
    }
  }
  message(sprintf("    %s 来源目录: %s", fs$id, ifelse(is.na(picked_dir), "未找到!", picked_dir)))
  for (ext in names(ext_dest_map)) {
    src_f <- if (!is.na(picked_dir)) file.path(picked_dir, paste0(fs$stem, ".", ext)) else NA_character_
    dst_f <- file.path(ext_dest_map[[ext]], paste0(fs$id, "_", fs$stem, ".", ext))
    copy_safe(src_f, dst_f, paste0("Main_Figure_", toupper(ext)))
  }
}

# ------------------------------------------------------------------------------
# 4. 生成机读基因集附件 (GMT + NHAS CSV) 与完整附表集合 (CSV + 合并 Excel)
# ------------------------------------------------------------------------------
message(">>> [4/6] 正在生成机读基因集附件与全部 Supplementary Tables...")

sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")
tcga_obj <- readRDS("01_clean_data/tcga_lihc_clean.rds")
h18_obj  <- readRDS("01_clean_data/hccdb18_clean.rds")
gse_obj  <- readRDS("01_clean_data/gse14520_clean.rds")

# 4.1 导出标准 GMT 与 NHAS 22 基因清单
gmt_orig_path <- file.path(dir_sup, "Supplementary_Data_S1_Original_24_GeneSets.gmt")
gmt_s3_path   <- file.path(dir_sup, "Supplementary_Data_S2_S3_Purified_GeneSets.gmt")
nhas_csv_path <- file.path(dir_sup, "Supplementary_Data_S3_NHAS_22_Genes.csv")

writeLines(
  vapply(names(sigs_obj$sig_list), function(nm) {
    paste(c(nm, "Original_24_HCC_Admixture_Benchmark", sigs_obj$sig_list[[nm]]), collapse = "\t")
  }, character(1)),
  gmt_orig_path
)
writeLines(
  vapply(names(res5_obj$purified_sig_list), function(nm) {
    paste(c(nm, "S3_scRNA_Purified_or_Unchanged", res5_obj$purified_sig_list[[nm]]), collapse = "\t")
  }, character(1)),
  gmt_s3_path
)
nhas_df <- tibble(Rank = seq_along(res3_obj$nhas_genes), Gene_Symbol = res3_obj$nhas_genes)
write.csv(nhas_df, nhas_csv_path, row.names = FALSE)

log_item("Supplementary_Data", gmt_orig_path, "Generated from immune_signatures.rds", "GENERATED")
log_item("Supplementary_Data", gmt_s3_path,   "Generated from 05_correction_benchmark_results.rds", "GENERATED")
log_item("Supplementary_Data", nhas_csv_path, "Generated from 03_scrna_mechanism_and_nhas.rds", "GENERATED")

# 4.2 构建各张 Supplementary Tables (统一收集到 list 用于写 CSV 及 Excel)
supp_tables <- list()

# S1: 基因集定义、S3 删减统计及三平台基因覆盖率
plat_mats <- list(TCGA_LIHC = tcga_obj$T_tpm, ICGC_LIRI_JP = h18_obj$T_tpm, GSE14520 = gse_obj$T_lin)
cov_wide <- map_dfr(names(sigs_obj$sig_list), function(sn) {
  g_orig <- sigs_obj$sig_list[[sn]]
  g_s3   <- if (sn %in% names(res5_obj$purified_sig_list)) res5_obj$purified_sig_list[[sn]] else g_orig
  tibble(
    Signature       = sn,
    Orig_Genes      = length(g_orig),
    S3_Kept_Genes   = length(g_s3),
    TCGA_Detected   = sum(g_orig %in% rownames(plat_mats$TCGA_LIHC)),
    ICGC_Detected   = sum(g_orig %in% rownames(plat_mats$ICGC_LIRI_JP)),
    GSE14520_Detect = sum(g_orig %in% rownames(plat_mats$GSE14520)),
    GSE14520_CovPct = round(sum(g_orig %in% rownames(plat_mats$GSE14520)) / length(g_orig) * 100, 1),
    Original_Genes_List = paste(g_orig, collapse = ";"),
    S3_Purified_List    = paste(g_s3, collapse = ";")
  )
}) %>%
  left_join(sigs_obj$sig_meta %>% select(Signature, Axis), by = "Signature") %>%
  relocate(Axis, .after = Signature)

supp_tables[["Table_S1_GeneSets_Coverage"]] <- cov_wide
supp_tables[["Table_S2_NHAS_22_Genes"]]     <- nhas_df

# 读取已有分析结果 CSV（纯只读加载）
csv_sources <- list(
  "Table_S3_Full_Atlas_24Sigs"      = "03_results/Table1_HCC_Immune_Admixture_Atlas.csv",
  "Table_S4_scRNA_TriArm_Reweight"  = "03_results/03b_scrna_physiological_triarm.csv",
  "Table_S5_Spatial_All7_QC"        = "03_results/revision_20261123/07b_spatial_all7_slides_qc.csv",
  "Table_S6_Spatial_All7_Summary"   = "03_results/revision_20261123/07b_spatial_all7_sig_summary.csv",
  "Table_S7_Bootstrap_95CI"         = "03_results/revision_20261123/07a_patient_bootstrap_95CI.csv",
  "Table_S8_MatchedNull_Grand"      = "03_results/revision_20261123/07a_random_deletion_null_grand.csv",
  "Table_S9_MatchedNull_BySig"      = "03_results/revision_20261123/07a_random_deletion_null_bysig.csv",
  "Table_S10_External_S2_Transfer"  = "03_results/revision_20261123/07a_external_S2_validation.csv",
  "Table_S11_Interaction_Quartile"  = "03_results/revision_20261123/07b_interaction_sign_diagnosis.csv",
  "Table_S12_Stromal_Quartiles"     = "03_results/revision_20261123/07b_stromal_quartile_HRs.csv",
  "Table_S13_Clinical_Audit"        = "03_results/revision_20261123/07b_clinical_covariate_audit.csv",
  "Table_S14_Multivar_Meta_Cox"     = "03_results/revision_20261123/07b_multivar_pruned_meta.csv",
  "Table_S15_Multivar_Cohort_Cox"   = "03_results/revision_20261123/07b_multivar_pruned_cohort.csv"
)

for (tname in names(csv_sources)) {
  fp <- csv_sources[[tname]]
  if (file.exists(fp)) {
    supp_tables[[tname]] <- read.csv(fp, check.names = FALSE)
  } else {
    log_item("Supplementary_Table", file.path(dir_sup_csv, paste0(tname, ".csv")), fp, "MISSING")
  }
}

# 将全部附表写入 Tables_CSV/
for (tname in names(supp_tables)) {
  out_csv <- file.path(dir_sup_csv, paste0(tname, ".csv"))
  write.csv(supp_tables[[tname]], out_csv, row.names = FALSE)
  src_lbl <- if (tname %in% names(csv_sources)) csv_sources[[tname]] else "Generated from cleaned RDS"
  log_item("Supplementary_Table_CSV", out_csv, src_lbl, "GENERATED")
}

# 自动尝试打包为单份美观多工作表 Excel (Supplementary_Tables.xlsx)
xlsx_path <- file.path(dir_sup, "Supplementary_Tables.xlsx")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(
    fontName = "Times New Roman", fontSize = 10, fontColour = "#FFFFFF",
    fgFill = "#1B365D", textDecoration = "bold", halign = "center", valign = "center",
    border = "TopBottom"
  )
  body_style <- openxlsx::createStyle(fontName = "Times New Roman", fontSize = 9.5)
  
  for (sn in names(supp_tables)) {
    sh <- substr(sn, 1, 31) # Excel Sheet 名长度上限 31 字符
    openxlsx::addWorksheet(wb, sh)
    dat <- supp_tables[[sn]]
    openxlsx::writeData(wb, sh, dat, headerStyle = hdr_style)
    openxlsx::addStyle(wb, sh, body_style, rows = 2:(nrow(dat) + 1), cols = 1:ncol(dat), gridExpand = TRUE)
    openxlsx::freezePane(wb, sh, firstRow = TRUE)
    openxlsx::setColWidths(wb, sh, cols = 1:ncol(dat), widths = "auto")
  }
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  log_item("Supplementary_Excel", xlsx_path, "Compiled via openxlsx", "GENERATED")
  message("    已生成多工作表附表: ", xlsx_path)
} else if (requireNamespace("writexl", quietly = TRUE)) {
  names(supp_tables) <- substr(names(supp_tables), 1, 31)
  writexl::write_xlsx(supp_tables, xlsx_path)
  log_item("Supplementary_Excel", xlsx_path, "Compiled via writexl", "GENERATED")
  message("    已生成多工作表附表: ", xlsx_path)
} else {
  message("    [提示] 当前环境未安装 openxlsx/writexl，已将全部 ", length(supp_tables),
          " 张附表以 CSV 保存至 ", dir_sup_csv, "（如需合并 Excel，可运行 install.packages('writexl') 后重跑）。")
}

# ------------------------------------------------------------------------------
# 5. 复制核心可复现代码与运行环境记录 (只复制，不改动 02_scripts/)
# ------------------------------------------------------------------------------
message(">>> [5/6] 正在复制可复现 R 脚本与环境版本信息...")

all_scripts <- sort(list.files("02_scripts", pattern = "\\.[Rr]$", full.names = TRUE))
for (sc in all_scripts) {
  copy_safe(sc, file.path(dir_code, basename(sc)), "Reproducibility_Code")
}

# 保存 R 环境快照
sess_file <- file.path(dir_repro, "R_sessionInfo.txt")
pkg_file  <- file.path(dir_repro, "R_package_versions.csv")
writeLines(capture.output(sessionInfo()), sess_file)

key_pkgs <- c("GSVA", "Seurat", "survival", "IOBR", "msigdbr", "tidyverse", "patchwork", "Matrix", "data.table")
pkg_tbl <- tibble(
  Package = key_pkgs,
  Version = vapply(key_pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "Not_Installed"
  }, character(1))
)
write.csv(pkg_tbl, pkg_file, row.names = FALSE)
log_item("Reproducibility_Env", sess_file, "sessionInfo()", "GENERATED")
log_item("Reproducibility_Env", pkg_file,  "packageVersion()", "GENERATED")

# ------------------------------------------------------------------------------
# 6. 生成投稿清单 Manifest + MD5 校验码 + 一键压缩 ZIP 备份
# ------------------------------------------------------------------------------
message(">>> [6/6] 正在生成清单与 MD5 校验码，并打包同级 ZIP 归档...")

manifest_df <- bind_rows(manifest_rows)
all_pkg_files <- sort(list.files(pkg_root, recursive = TRUE, full.names = TRUE))
all_pkg_files <- all_pkg_files[!dir.exists(all_pkg_files)]
md5_vec <- tools::md5sum(all_pkg_files)

manifest_txt <- file.path(pkg_root, "00_Submission_Manifest_and_MD5.txt")
lines_out <- c(
  "==============================================================================",
  "JOURNAL SUBMISSION PACKAGE MANIFEST",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Package Directory: ", normalizePath(pkg_root)),
  "Note: All files were strictly COPIED or freshly generated; original project",
  "      directories were NOT modified or moved.",
  "==============================================================================",
  "",
  "--- [1] SUMMARY BY STATUS ---",
  capture.output(print(as.data.frame(manifest_df %>% count(Category, Status)), row.names = FALSE)),
  "",
  "--- [2] DETAILED FILE PROVENANCE ---",
  capture.output(print(as.data.frame(manifest_df), row.names = FALSE)),
  "",
  "--- [3] MD5 CHECKSUMS ---",
  paste(unname(md5_vec), sub(paste0("^", pkg_root, "/"), "", names(md5_vec)), sep = "  ")
)
writeLines(lines_out, manifest_txt)

# 生成同级 ZIP 压缩包（不删除文件夹，方便备份或直接传输）
zip_out <- paste0(pkg_root, ".zip")
zip_ok <- tryCatch({
  if (file.exists(zip_out)) file.remove(zip_out)
  status <- system2("zip", args = c("-r", "-q", zip_out, pkg_root))
  status == 0
}, error = function(e) FALSE)

message("\n==============================================================================")
message(">>> 投稿包全部打包完成！原有项目目录 100% 保持原样未改动。")
message("    新建投稿文件夹: ", normalizePath(pkg_root))
if (zip_ok && file.exists(zip_out)) {
  message(sprintf("    同级 ZIP 归档 : %s (%.1f MB)", normalizePath(zip_out), file.size(zip_out) / 1024^2))
}
message("==============================================================================")
print(as.data.frame(manifest_df %>% count(Category, Status)), row.names = FALSE)

if (any(manifest_df$Status %in% c("MISSING", "ACTION_NEEDED"))) {
  message("\n[请留意以下待确认项]:")
  print(as.data.frame(manifest_df %>% filter(Status %in% c("MISSING", "ACTION_NEEDED"))), row.names = FALSE)
}
