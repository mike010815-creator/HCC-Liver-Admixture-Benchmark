# ==============================================================================
# 脚本名称: 01_clean_tcga.R
# 核心任务:
#   1. 读取 TCGA-LIHC STAR-TPM，还原为线性 TPM 并按样本重标化至 1e6
#   2. Ensembl ID 转标准 Gene Symbol 并按最高均值去重
#   3. 精准提取配对癌 (01) 与癌旁 (11)，以及全队列原发癌 + 生存信息
#   4. 组装 4 大生物学轴的 HCC 免疫转录组特征库 (MSigDB + IOBR + 经典文献)
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(msigdbr)
  library(IOBR)
})

# ------------------------------------------------------------------------------
# 1. 检查原始文件是否存在
# ------------------------------------------------------------------------------
f_tpm   <- "00_raw_data/TCGA_LIHC/TCGA-LIHC.star_tpm.tsv"
f_probe <- "00_raw_data/TCGA_LIHC/gencode.v36.annotation.gtf.gene.probemap"
f_surv  <- "00_raw_data/TCGA_LIHC/TCGA-LIHC.survival.tsv"
f_clin  <- "00_raw_data/TCGA_LIHC/TCGA-LIHC.clinical.tsv"

stopifnot(
  "找不到 star_tpm.tsv，请检查工作目录和文件名" = file.exists(f_tpm),
  "找不到 probemap 文件" = file.exists(f_probe),
  "找不到 survival.tsv 文件" = file.exists(f_surv),
  "找不到 clinical.tsv 文件" = file.exists(f_clin)
)
message(">>> [1/6] 原始文件路径检查全部通过！")

# ------------------------------------------------------------------------------
# 2. 读取 TPM 矩阵并还原回线性尺度
# ------------------------------------------------------------------------------
tpm_raw  <- fread(f_tpm, data.table = FALSE)
probemap <- fread(f_probe, data.table = FALSE)

# 第 1 列为 Ensembl_ID，其余列为样本表达值
ens_ids  <- tpm_raw[[1]]
expr_mat <- as.matrix(tpm_raw[, -1, drop = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- ens_ids

max_val <- max(expr_mat, na.rm = TRUE)
message(sprintf(">>> [2/6] 原始矩阵维度: %d 基因 x %d 样本 | 最大值: %.2f",
                nrow(expr_mat), ncol(expr_mat), max_val))

# UCSC Xena GDC STAR-TPM 默认单位为 log2(TPM + 1)，最大值通常在 15-20 之间
if (max_val < 30) {
  message("    检测到 log2(TPM+1) 尺度，正在还原为线性 TPM...")
  tpm_linear <- 2^expr_mat - 1
  tpm_linear[tpm_linear < 0] <- 0
} else {
  message("    数据已是线性尺度，直接使用。")
  tpm_linear <- expr_mat
}

# ------------------------------------------------------------------------------
# 3. Ensembl ID -> Gene Symbol 转换与去重
# ------------------------------------------------------------------------------
# 优先精确匹配 Ensembl ID；若匹配极少则去掉小数点版本号再匹配
idx <- match(rownames(tpm_linear), probemap$id)
if (sum(!is.na(idx)) < 1000) {
  idx <- match(sub("\\..*$", "", rownames(tpm_linear)),
               sub("\\..*$", "", probemap$id))
}

gene_sym <- probemap$gene[idx]
keep_valid <- !is.na(gene_sym) & gene_sym != "" & !grepl("^\\s*$", gene_sym)

tpm_valid <- tpm_linear[keep_valid, , drop = FALSE]
sym_valid <- gene_sym[keep_valid]

# 重复基因名保留全队列平均表达量最高的那一行
row_means <- rowMeans(tpm_valid, na.rm = TRUE)
ord <- order(sym_valid, -row_means)
tpm_valid <- tpm_valid[ord, , drop = FALSE]
sym_valid <- sym_valid[ord]

keep_unique <- !duplicated(sym_valid)
tpm_dedup <- tpm_valid[keep_unique, , drop = FALSE]
rownames(tpm_dedup) <- sym_valid[keep_unique]

# 过滤掉在所有样本中表达量恒为 0 的基因
tpm_dedup <- tpm_dedup[rowSums(tpm_dedup) > 0, , drop = FALSE]

# 检查去重后每列总和，并重标化为严格 1e6（确保后续线性混合时总转录本质量严格守恒）
raw_colsums <- summary(colSums(tpm_dedup))
message(sprintf(">>> [3/6] 基因转换去重完成: 保留 %d 个唯一基因 | 重标化前样本总 TPM 中位数: %.0f",
                nrow(tpm_dedup), raw_colsums["Median"]))

tpm_norm <- sweep(tpm_dedup, 2, colSums(tpm_dedup), FUN = "/") * 1e6

# ------------------------------------------------------------------------------
# 4. 提取配对癌 (01) 与癌旁 (11)，以及全队列原发癌 + 生存信息
# ------------------------------------------------------------------------------
sample_df <- tibble(
  sample_id  = colnames(tpm_norm),
  patient_id = substr(colnames(tpm_norm), 1, 12),
  type_code  = substr(colnames(tpm_norm), 14, 15)
) %>%
  arrange(sample_id)  # 排序确保同一患者若同时有 01A 和 01B，稳定优先取 01A

tumor_df  <- sample_df %>% filter(type_code == "01") %>% distinct(patient_id, .keep_all = TRUE)
normal_df <- sample_df %>% filter(type_code == "11") %>% distinct(patient_id, .keep_all = TRUE)

# 配对交集
paired_patients <- sort(intersect(tumor_df$patient_id, normal_df$patient_id))
paired_T_df <- tumor_df  %>% filter(patient_id %in% paired_patients) %>% arrange(patient_id)
paired_N_df <- normal_df %>% filter(patient_id %in% paired_patients) %>% arrange(patient_id)

# 强制断言：两张表的患者顺序必须 100% 一一对应
stopifnot(identical(paired_T_df$patient_id, paired_N_df$patient_id))

T_tpm <- tpm_norm[, paired_T_df$sample_id, drop = FALSE]
N_tpm <- tpm_norm[, paired_N_df$sample_id, drop = FALSE]

# 全队列原发肿瘤（用于第 4-5 关真实生存分析）
all_T_tpm <- tpm_norm[, tumor_df$sample_id, drop = FALSE]

# 对齐临床与生存信息表
surv_df <- fread(f_surv, data.table = FALSE)
clin_df <- fread(f_clin, data.table = FALSE)
colnames(surv_df)[1] <- "sample_id"
colnames(clin_df)[1] <- "sample_id"

all_T_meta <- tumor_df %>%
  left_join(surv_df, by = "sample_id") %>%
  left_join(clin_df, by = "sample_id")

n_valid_os <- sum(!is.na(all_T_meta$OS) & !is.na(all_T_meta$OS.time) & all_T_meta$OS.time > 0)

message(sprintf(">>> [4/6] 样本提取完成: 配对癌-癌旁 = %d 对 | 全队列原发癌 = %d 例 (有效 OS = %d 例)",
                ncol(T_tpm), ncol(all_T_tpm), n_valid_os))

# ------------------------------------------------------------------------------
# 5. 生物学常识自检（确认癌与癌旁绝对没有贴反标签）
# ------------------------------------------------------------------------------
qc_genes <- c("ALB", "APOA2", "CYP3A4", "GPC3", "AFP", "MKI67")
qc_genes <- intersect(qc_genes, rownames(T_tpm))

qc_table <- tibble(
  Gene = qc_genes,
  Expected_Higher_In = ifelse(qc_genes %in% c("ALB", "APOA2", "CYP3A4"), "Normal (癌旁)", "Tumor (肝癌)"),
  Tumor_Median_TPM   = round(apply(T_tpm[qc_genes, , drop = FALSE], 1, median), 1),
  Normal_Median_TPM  = round(apply(N_tpm[qc_genes, , drop = FALSE], 1, median), 1)
)
message(">>> [5/6] 生物学质控自检表（请核对中位数方向是否与预期一致）:")
print(as.data.frame(qc_table), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. 组装 4 大生物学轴的 HCC 免疫转录组特征库
# ------------------------------------------------------------------------------
# 6.1 从 MSigDB 提取权威 Hallmark 免疫/炎症/基质排斥通路（兼容新旧版 msigdbr）
msig_h <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
get_hallmark <- function(gs_name_str) {
  unique(msig_h$gene_symbol[msig_h$gs_name == gs_name_str])
}

# 6.2 从 IOBR 安全提取已发表 TME 特征（自动不区分大小写匹配，防止不同版本报错）
safe_load_iobr <- function(name) {
  res <- tryCatch({
    if (exists("load_data", envir = asNamespace("IOBR"), mode = "function")) {
      IOBR::load_data(name)
    } else NULL
  }, error = function(e) NULL)
  
  if (is.null(res)) {
    res <- tryCatch({
      e <- new.env()
      data(list = name, package = "IOBR", envir = e)
      e[[name]]
    }, error = function(e) NULL)
  }
  return(res)
}

signature_tme        <- safe_load_iobr("signature_tme")
signature_collection <- safe_load_iobr("signature_collection")
iobr_pool            <- c(signature_tme, signature_collection)
iobr_names_lower <- tolower(names(iobr_pool))

fetch_iobr <- function(candidates) {
  for (cand in candidates) {
    hit <- which(iobr_names_lower == tolower(cand))
    if (length(hit) > 0) return(unique(iobr_pool[[hit[1]]]))
  }
  warning("IOBR 中未匹配到特征: ", paste(candidates, collapse = " / "))
  return(NULL)
}

# 6.3 汇总 4 大轴共 22 个核心特征
sig_raw <- list(
  # --- 轴 1：抗肿瘤适应性免疫 (Adaptive Anti-tumor) ---
  CYT_Rooney             = c("GZMA", "PRF1"),                                  # Rooney et al., Cell 2015
  T_cell_inflamed_GEP    = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1",
                             "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E",
                             "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10",
                             "STAT1", "TIGIT"),                                # Ayers et al., JCI 2017
  CD8_T_effector         = fetch_iobr(c("CD_8_T_effector")),                   # Mariathasan et al., Nature 2018
  HM_IFN_gamma           = get_hallmark("HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  HM_IFN_alpha           = get_hallmark("HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  HM_Allograft_Rejection = get_hallmark("HALLMARK_ALLOGRAFT_REJECTION"),
  
  # --- 轴 2：免疫抑制、耗竭与髓系微环境 (Suppression / Myeloid) ---
  Immune_Checkpoint      = fetch_iobr(c("Immune_Checkpoint")),
  T_cell_exhaustion      = fetch_iobr(c("T_cell_exhaustion", "Exhausted_CD8_T_cell", "T_cell_co_inhibition")),
  Treg_Signature         = fetch_iobr(c("Tregs_Bindea_et_al", "Treg_Peng_et_al", "T_regulatory_cells")),
  MDSC_Peng              = fetch_iobr(c("MDSC_Peng_et_al", "MDSC")),
  TAM_Peng               = fetch_iobr(c("TAM_Peng_et_al", "Macrophages_M2_CIBERSORT", "M2_macrophage")),
  HM_IL2_STAT5           = get_hallmark("HALLMARK_IL2_STAT5_SIGNALING"),
  
  # --- 轴 3：固有炎症与肝微环境易串扰轴 (Innate Inflammation / Complement) ---
  HM_Complement          = get_hallmark("HALLMARK_COMPLEMENT"),
  HM_Inflammatory_Resp   = get_hallmark("HALLMARK_INFLAMMATORY_RESPONSE"),
  HM_IL6_JAK_STAT3       = get_hallmark("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  HM_TNFA_via_NFKB       = get_hallmark("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  
  # --- 轴 4：基质纤维化与免疫排斥轴 (Stromal / Exclusion) ---
  Pan_F_TBRs             = fetch_iobr(c("Pan_F_TBRs")),                        # Mariathasan et al., Nature 2018
  CAF_Peng               = fetch_iobr(c("CAF_Peng_et_al", "CAF")),
  WNT_target             = fetch_iobr(c("WNT_target")),
  HM_TGF_beta            = get_hallmark("HALLMARK_TGF_BETA_SIGNALING"),
  HM_EMT                 = get_hallmark("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  HM_Wnt_Beta_Catenin    = get_hallmark("HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
  HM_Angiogenesis        = get_hallmark("HALLMARK_ANGIOGENESIS")
)

# 去掉未能匹配到的空条目（如果有）
sig_list <- sig_raw[!sapply(sig_raw, is.null)]

# 为每个特征标注所属生物学轴并统计在 TCGA 矩阵中的基因检出率
axis_map <- c(
  CYT_Rooney = "1_Adaptive_Immunity", T_cell_inflamed_GEP = "1_Adaptive_Immunity",
  CD8_T_effector = "1_Adaptive_Immunity", HM_IFN_gamma = "1_Adaptive_Immunity",
  HM_IFN_alpha = "1_Adaptive_Immunity", HM_Allograft_Rejection = "1_Adaptive_Immunity",
  Immune_Checkpoint = "2_Suppression_Myeloid", T_cell_exhaustion = "2_Suppression_Myeloid",
  Treg_Signature = "2_Suppression_Myeloid", MDSC_Peng = "2_Suppression_Myeloid",
  TAM_Peng = "2_Suppression_Myeloid", HM_IL2_STAT5 = "2_Suppression_Myeloid",
  HM_Complement = "3_Innate_Inflammation", HM_Inflammatory_Resp = "3_Innate_Inflammation",
  HM_IL6_JAK_STAT3 = "3_Innate_Inflammation", HM_TNFA_via_NFKB = "3_Innate_Inflammation",
  Pan_F_TBRs = "4_Stromal_Exclusion", CAF_Peng = "4_Stromal_Exclusion",
  WNT_target = "4_Stromal_Exclusion", HM_TGF_beta = "4_Stromal_Exclusion",
  HM_EMT = "4_Stromal_Exclusion", HM_Wnt_Beta_Catenin = "4_Stromal_Exclusion",
  HM_Angiogenesis = "4_Stromal_Exclusion"
)

sig_meta <- tibble(
  Signature      = names(sig_list),
  Axis           = unname(axis_map[names(sig_list)]),
  Total_Genes    = sapply(sig_list, length),
  Detected_TCGA  = sapply(sig_list, function(g) sum(g %in% rownames(T_tpm))),
  Detect_Rate    = round(Detected_TCGA / Total_Genes * 100, 1)
)

message(">>> [6/6] 免疫特征库构建完成，检出率汇总如下:")
print(as.data.frame(sig_meta), row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. 保存全部清洗结果至 01_clean_data/ 与 03_results/
# ------------------------------------------------------------------------------
saveRDS(
  list(
    T_tpm       = T_tpm,
    N_tpm       = N_tpm,
    paired_meta = tibble(patient_id = paired_patients,
                         T_sample = paired_T_df$sample_id,
                         N_sample = paired_N_df$sample_id),
    all_T_tpm   = all_T_tpm,
    all_T_meta  = all_T_meta
  ),
  file = "01_clean_data/tcga_lihc_clean.rds"
)

saveRDS(
  list(sig_list = sig_list, sig_meta = sig_meta),
  file = "01_clean_data/immune_signatures.rds"
)

write.csv(sig_meta, "03_results/01_signature_library_summary.csv", row.names = FALSE)
message(">>> 全部完成！已生成 01_clean_data/tcga_lihc_clean.rds 与 immune_signatures.rds。")
