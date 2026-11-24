# ==============================================================================
# 脚本名称: 02d_external_bulk_val.R
# 核心任务:
#   1. 修复单基因 Treg_Signature，统一 24 个特征库，并补跑 TCGA-LIHC 24 特征基准
#   2. 清洗并配对提取 ICGC-LIRI-JP (HCCDB18, RNA-seq) -> 保存 hccdb18_clean.rds
#   3. 清洗并配对提取 GSE14520 (GPL3921, Microarray)  -> 保存 gse14520_clean.rds
#   4. 在三大独立队列上运行统一的 0%–50% 线性混合模拟，评估跨队列/跨平台一致性
#   5. 导出三队列一致性散点图与热图 (Fig2C_MultiCohort_Consistency)
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GSVA)
  library(IOBR)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 修复 Treg_Signature (确保基因数 >= 3)，统一 24 个特征库
# ------------------------------------------------------------------------------
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
sig_list <- sigs_obj$sig_list
sig_meta <- sigs_obj$sig_meta

iobr_pool <- c(IOBR::load_data("signature_tme"), IOBR::load_data("signature_collection"))
treg_cands <- c("CD4_c1_Treg", "Treg_Rooney_et_al", "Treg_Danaher_et_al", "Treg_Bindea_et_al")
treg_lens  <- sapply(treg_cands, function(k) if (k %in% names(iobr_pool)) length(unique(iobr_pool[[k]])) else 0)
message(">>> [1/6] 本机 IOBR Treg 候选特征基因数: ",
        paste(names(treg_lens), treg_lens, sep = "=", collapse = " | "))

# 优先选用基因数 >= 3 且最具肝癌/泛癌代表性的多基因 Treg 特征
valid_treg <- names(treg_lens)[treg_lens >= 3]
if (length(valid_treg) > 0) {
  chosen_treg <- valid_treg[1]
  sig_list$Treg_Signature <- unique(iobr_pool[[chosen_treg]])
  message(sprintf("    已将 Treg_Signature 替换为 %s (%d 个基因)。",
                  chosen_treg, length(sig_list$Treg_Signature)))
} else {
  # 兜底: 若 IOBR 所有 Treg 均 <3 基因，则采用标准多基因经典 Treg 标记集合
  sig_list$Treg_Signature <- c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF4", "TNFRSF18", "IKZF2", "CCR8")
  message("    IOBR 候选均 <3 基因，已替换为标准 8 基因 Treg 集合。")
}

sig_meta$Total_Genes <- sapply(sig_meta$Signature, function(s) length(sig_list[[s]]))
saveRDS(list(sig_list = sig_list, sig_meta = sig_meta),
        file = "01_clean_data/immune_signatures.rds")

# ------------------------------------------------------------------------------
# 2. 清洗外部队列 1: ICGC-LIRI-JP (HCCDB18)
# ------------------------------------------------------------------------------
message("\n>>> [2/6] 正在清洗 ICGC-LIRI-JP (HCCDB18)...")
h18_exp_raw  <- fread("00_raw_data/HCCDB18_ICGC/HCCDB18_mRNA_level3.txt", data.table = FALSE)
h18_samp_raw <- fread("00_raw_data/HCCDB18_ICGC/HCCDB18.sample.txt", data.table = FALSE)
h18_pat_raw  <- fread("00_raw_data/HCCDB18_ICGC/HCCDB18.patient.txt", data.table = FALSE)

# 转置 sample.txt 与 patient.txt（HCCDB 为“行=属性，列=样本/患者”）
transpose_hccdb <- function(df, id_colname) {
  attr_names <- make.unique(as.character(df[[1]]))
  mat <- t(as.matrix(df[, -1, drop = FALSE]))
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(out) <- attr_names
  out[[id_colname]] <- rownames(out)
  rownames(out) <- NULL
  as_tibble(out)
}

h18_samp <- transpose_hccdb(h18_samp_raw, "sample_id")
h18_pat  <- transpose_hccdb(h18_pat_raw,  "patient_id")

message("    HCCDB18.sample 列名: ", paste(colnames(h18_samp), collapse = ", "))
message("    HCCDB18.patient 列名: ", paste(colnames(h18_pat), collapse = ", "))

# 构建基因表达矩阵并去重
sym_col <- grep("symbol", colnames(h18_exp_raw), ignore.case = TRUE)[1]
h18_syms <- as.character(h18_exp_raw[[sym_col]])
h18_mat  <- as.matrix(h18_exp_raw[, (sym_col + 1):ncol(h18_exp_raw), drop = FALSE])
mode(h18_mat) <- "numeric"

# 过滤空基因名与含 NA 过多的行，对少量缺失填 0（或仅保留完整基因）
keep_g <- !is.na(h18_syms) & h18_syms != "" & rowSums(is.na(h18_mat)) == 0
h18_mat  <- h18_mat[keep_g, , drop = FALSE]
h18_syms <- h18_syms[keep_g]

ord <- order(h18_syms, -rowMeans(h18_mat))
h18_mat  <- h18_mat[ord, , drop = FALSE]
h18_syms <- h18_syms[ord]
h18_mat  <- h18_mat[!duplicated(h18_syms), , drop = FALSE]
rownames(h18_mat) <- h18_syms[!duplicated(h18_syms)]

# 检查数值尺度并还原线性尺度 + 每列重标化至 1e6
h18_min <- min(h18_mat); h18_max <- max(h18_mat)
message(sprintf("    HCCDB18 原始数值范围: [%.2f, %.2f]", h18_min, h18_max))
if (h18_max < 30) {
  h18_lin <- if (h18_min >= 0) pmax(2^h18_mat - 1, 0) else 2^h18_mat
} else {
  h18_lin <- pmax(h18_mat, 0)
}
h18_lin <- h18_lin[rowSums(h18_lin) > 0, , drop = FALSE]
h18_lin <- sweep(h18_lin, 2, colSums(h18_lin), "/") * 1e6

# 配对提取 HCC 与 Adjacent
pat_link_col <- grep("PATIENT", colnames(h18_samp), ignore.case = TRUE, value = TRUE)[1]
type_col     <- grep("^TYPE$",  colnames(h18_samp), ignore.case = TRUE, value = TRUE)[1]

h18_samp_in <- h18_samp %>%
  filter(sample_id %in% colnames(h18_lin)) %>%
  rename(pat_key = all_of(pat_link_col), tissue_type = all_of(type_col))

h18_T_df <- h18_samp_in %>% filter(grepl("^HCC$|Tumor", tissue_type, ignore.case = TRUE)) %>%
  arrange(sample_id) %>% distinct(pat_key, .keep_all = TRUE)
h18_N_df <- h18_samp_in %>% filter(grepl("Adj|Normal|Non", tissue_type, ignore.case = TRUE)) %>%
  arrange(sample_id) %>% distinct(pat_key, .keep_all = TRUE)

h18_pairs <- sort(intersect(h18_T_df$pat_key, h18_N_df$pat_key))
h18_T_pair <- h18_T_df %>% filter(pat_key %in% h18_pairs) %>% arrange(pat_key)
h18_N_pair <- h18_N_df %>% filter(pat_key %in% h18_pairs) %>% arrange(pat_key)
stopifnot(identical(h18_T_pair$pat_key, h18_N_pair$pat_key))

h18_T_tpm <- h18_lin[, h18_T_pair$sample_id, drop = FALSE]
h18_N_tpm <- h18_lin[, h18_N_pair$sample_id, drop = FALSE]
h18_all_T <- h18_lin[, h18_T_df$sample_id,   drop = FALSE]

message(sprintf("    HCCDB18 清洗完成: %d 基因 | 配对 = %d 对 | 全原发癌 = %d 例 | ALB 中位(T vs N): %.0f vs %.0f",
                nrow(h18_T_tpm), ncol(h18_T_tpm), ncol(h18_all_T),
                median(h18_T_tpm["ALB", ]), median(h18_N_tpm["ALB", ])))

saveRDS(
  list(T_tpm = h18_T_tpm, N_tpm = h18_N_tpm,
       all_T_tpm = h18_all_T, samp_meta = h18_samp_in, pat_meta = h18_pat),
  file = "01_clean_data/hccdb18_clean.rds"
)

# ------------------------------------------------------------------------------
# 3. 清洗外部队列 2: GSE14520 (GPL3921 芯片队列)
# ------------------------------------------------------------------------------
message("\n>>> [3/6] 正在清洗 GSE14520 (GPL3921)...")
f_gse_mat <- "00_raw_data/GSE14520/GSE14520-GPL3921_series_matrix.txt.gz"
f_gse_sup <- "00_raw_data/GSE14520/GSE14520_Extra_Supplement.txt.gz"
f_gse_gpl <- "00_raw_data/GSE14520/GPL3921_annot.txt"

# 精确定位 series_matrix 中 ID_REF 表头行号，避开顶部 ! 注释
con <- gzfile(f_gse_mat, "rt")
head_lines <- readLines(con, n = 200)
close(con)
skip_mat <- which(grepl('^"?ID_REF"?\t', head_lines))[1] - 1
stopifnot("无法定位 GSE14520 的 ID_REF 表头" = !is.na(skip_mat))

gse_raw <- fread(f_gse_mat, skip = skip_mat, data.table = FALSE)
gse_raw <- gse_raw[!grepl("^!", gse_raw[[1]]) & gse_raw[[1]] != "", , drop = FALSE]
probe_ids <- gsub('^"|"$', "", as.character(gse_raw[[1]]))
gse_mat   <- as.matrix(gse_raw[, -1, drop = FALSE])
mode(gse_mat) <- "numeric"
rownames(gse_mat) <- probe_ids
colnames(gse_mat) <- gsub('^"|"$', "", colnames(gse_mat))

# 读取 GPL3921 探针注释表（自动定位 ID 表头并防引号解析报错）
gpl_head <- readLines(f_gse_gpl, n = 100)
skip_gpl <- which(grepl('^"?ID"?\t', gpl_head))[1] - 1
gpl_df   <- fread(f_gse_gpl, skip = skip_gpl, sep = "\t", quote = "", fill = TRUE, data.table = FALSE)

gpl_id_col  <- colnames(gpl_df)[1]
gpl_sym_col <- grep("^Gene[ _.]?Symbol$", colnames(gpl_df), ignore.case = TRUE, value = TRUE)[1]
stopifnot("GPL3921 中找不到 Gene Symbol 列" = !is.na(gpl_sym_col))

gpl_ids  <- gsub('^"|"$', "", as.character(gpl_df[[gpl_id_col]]))
gpl_syms <- gsub('^"|"$', "", as.character(gpl_df[[gpl_sym_col]]))

# 映射探针 -> Gene Symbol，剔除空值与多基因歧义探针 ("///")，重复基因保留高均值探针
p_idx  <- match(rownames(gse_mat), gpl_ids)
mapped <- trimws(gpl_syms[p_idx])
keep_p <- !is.na(mapped) & mapped != "" & mapped != "---" &
  !grepl("///", mapped) & rowSums(is.na(gse_mat)) == 0

gse_valid <- gse_mat[keep_p, , drop = FALSE]
sym_valid <- mapped[keep_p]

ord_g <- order(sym_valid, -rowMeans(gse_valid))
gse_valid <- gse_valid[ord_g, , drop = FALSE]
sym_valid <- sym_valid[ord_g]
gse_dedup <- gse_valid[!duplicated(sym_valid), , drop = FALSE]
rownames(gse_dedup) <- sym_valid[!duplicated(sym_valid)]

# Affymetrix RMA 为 log2 强度 -> 还原线性尺度 2^x 并每列重标化至 1e6
message(sprintf("    GSE14520 原始数值范围: [%.2f, %.2f]", min(gse_dedup), max(gse_dedup)))
gse_lin <- if (max(gse_dedup) < 30) 2^gse_dedup else pmax(gse_dedup, 0)
gse_lin <- sweep(gse_lin, 2, colSums(gse_lin), "/") * 1e6

# 关联 Extra_Supplement，精准提取配对 Tumor & Non-Tumor
gse_sup <- fread(f_gse_sup, data.table = FALSE)
message("    GSE14520_Extra_Supplement 列名: ", paste(colnames(gse_sup), collapse = ", "))

gsm_col  <- grep("Affy_GSM", colnames(gse_sup), ignore.case = TRUE, value = TRUE)[1]
tt_col   <- grep("Tissue",   colnames(gse_sup), ignore.case = TRUE, value = TRUE)[1]
lcs_col  <- grep("LCS",      colnames(gse_sup), ignore.case = TRUE, value = TRUE)[1]
pid_col  <- grep("^ID$",     colnames(gse_sup), ignore.case = TRUE, value = TRUE)[1]

sup_in <- gse_sup %>%
  filter(!is.na(.data[[gsm_col]]), .data[[gsm_col]] %in% colnames(gse_lin)) %>%
  mutate(
    sample_id   = .data[[gsm_col]],
    tissue_type = trimws(as.character(.data[[tt_col]])),
    # 同时尝试 ID 列与去掉末尾 A/B 的 LCS_ID 作为配对键，自动选配对数更多且唯一的键
    key_id  = trimws(as.character(.data[[pid_col]])),
    key_lcs = sub("[A-Za-z]$", "", trimws(as.character(.data[[lcs_col]])))
  )

t_rows <- sup_in %>% filter(grepl("^Tumor$", tissue_type, ignore.case = TRUE))
n_rows <- sup_in %>% filter(grepl("Non",     tissue_type, ignore.case = TRUE))

pairs_by_id  <- intersect(t_rows$key_id[t_rows$key_id != "" & !is.na(t_rows$key_id)],
                          n_rows$key_id[n_rows$key_id != "" & !is.na(n_rows$key_id)])
pairs_by_lcs <- intersect(t_rows$key_lcs[t_rows$key_lcs != "" & !is.na(t_rows$key_lcs)],
                          n_rows$key_lcs[n_rows$key_lcs != "" & !is.na(n_rows$key_lcs)])

use_col <- if (length(pairs_by_id) >= length(pairs_by_lcs)) "key_id" else "key_lcs"
common_keys <- sort(if (use_col == "key_id") pairs_by_id else pairs_by_lcs)

gse_T_pair <- t_rows %>% filter(.data[[use_col]] %in% common_keys) %>%
  distinct(.data[[use_col]], .keep_all = TRUE) %>% arrange(.data[[use_col]])
gse_N_pair <- n_rows %>% filter(.data[[use_col]] %in% common_keys) %>%
  distinct(.data[[use_col]], .keep_all = TRUE) %>% arrange(.data[[use_col]])
stopifnot(identical(gse_T_pair[[use_col]], gse_N_pair[[use_col]]))

gse_T_lin <- gse_lin[, gse_T_pair$sample_id, drop = FALSE]
gse_N_lin <- gse_lin[, gse_N_pair$sample_id, drop = FALSE]
gse_all_T <- gse_lin[, unique(t_rows$sample_id), drop = FALSE]

message(sprintf("    GSE14520 清洗完成: %d 基因 | 配对 (%s) = %d 对 | 全原发癌 = %d 例 | ALB 中位(T vs N): %.0f vs %.0f",
                nrow(gse_T_lin), use_col, ncol(gse_T_lin), ncol(gse_all_T),
                median(gse_T_lin["ALB", ]), median(gse_N_lin["ALB", ])))

saveRDS(
  list(T_lin = gse_T_lin, N_lin = gse_N_lin,
       all_T_lin = gse_all_T, meta = sup_in, pair_key = use_col),
  file = "01_clean_data/gse14520_clean.rds"
)

# ------------------------------------------------------------------------------
# 4. 统一函数: 在 3 大独立队列上运行 0%–50% 线性混合模拟 (ssGSEA + Zscore)
# ------------------------------------------------------------------------------
run_cohort_sim <- function(T_mat, N_mat, cohort_name, gsets) {
  alphas <- c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
  
  mix_logs <- lapply(alphas, function(a) {
    log2((1 - a) * T_mat + a * N_mat + 1)
  })
  
  # ssGSEA (纯单样本, normalize = FALSE)
  ss_list <- lapply(mix_logs, function(m) {
    p <- ssgseaParam(as.matrix(m), gsets, minSize = 2, normalize = FALSE)
    gsva(p, verbose = FALSE)
  })
  
  # 固定基线 Z-score
  m0 <- mix_logs[[1]]
  mu0 <- rowMeans(m0); sd0 <- apply(m0, 1, sd); sd0[sd0 == 0] <- NA
  zs_list <- lapply(mix_logs, function(m) {
    zm <- sweep(sweep(m, 1, mu0, "-"), 1, sd0, "/")
    out <- sapply(gsets, function(g) {
      gu <- intersect(g, rownames(zm))
      colMeans(zm[gu, , drop = FALSE], na.rm = TRUE)
    })
    t(out)
  })
  
  eval_list <- function(slist, method_lbl) {
    s0_mat <- slist[[1]]
    map_dfr(seq_along(alphas), function(k) {
      a <- alphas[k]; sa_mat <- slist[[k]]
      map_dfr(rownames(s0_mat), function(sig) {
        s0 <- as.numeric(s0_mat[sig, ]); sa <- as.numeric(sa_mat[sig, ])
        med0 <- median(s0)
        tibble(
          Cohort         = cohort_name,
          Method         = method_lbl,
          Signature      = sig,
          Alpha_Pct      = a * 100,
          Std_Drift      = mean(sa - s0) / sd(s0),
          Spearman_Rho   = cor(s0, sa, method = "spearman"),
          Flip_Fixed_Pct = mean((s0 >= med0) != (sa >= med0)) * 100
        )
      })
    })
  }
  
  bind_rows(eval_list(ss_list, "ssGSEA"), eval_list(zs_list, "Zscore"))
}

tcga <- readRDS("01_clean_data/tcga_lihc_clean.rds")

message("\n>>> [4/6] 正在三大队列上并行运行 0%–50% 混入模拟（约 20 秒）...")
sim_all <- bind_rows(
  run_cohort_sim(tcga$T_tpm, tcga$N_tpm, "TCGA_LIHC (n=50, RNA-seq)",      sig_list),
  run_cohort_sim(h18_T_tpm,  h18_N_tpm,  sprintf("ICGC_LIRI_JP (n=%d, RNA-seq)", ncol(h18_T_tpm)), sig_list),
  run_cohort_sim(gse_T_lin,  gse_N_lin,  sprintf("GSE14520 (n=%d, Microarray)",  ncol(gse_T_lin)), sig_list)
) %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature")

saveRDS(sim_all, "01_clean_data/02d_multicohort_sim_results.rds")
write.csv(sim_all, "03_results/02d_multicohort_admixture_metrics.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. 计算跨队列敏感性排序的一致性 (Spearman & Pearson r)
# ------------------------------------------------------------------------------
wide30_ss <- sim_all %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  mutate(C_short = case_when(
    grepl("TCGA", Cohort) ~ "TCGA",
    grepl("ICGC", Cohort) ~ "ICGC",
    TRUE                  ~ "GSE14520"
  )) %>%
  select(Signature, Axis, C_short, Std_Drift, Flip_Fixed_Pct) %>%
  pivot_wider(names_from = C_short, values_from = c(Std_Drift, Flip_Fixed_Pct)) %>%
  arrange(desc(Std_Drift_TCGA))

r_tcga_icgc_s <- cor(wide30_ss$Std_Drift_TCGA, wide30_ss$Std_Drift_ICGC, method = "spearman")
r_tcga_gse_s  <- cor(wide30_ss$Std_Drift_TCGA, wide30_ss$Std_Drift_GSE14520, method = "spearman")
r_icgc_gse_s  <- cor(wide30_ss$Std_Drift_ICGC, wide30_ss$Std_Drift_GSE14520, method = "spearman")

wide30_zs <- sim_all %>%
  filter(Method == "Zscore", Alpha_Pct == 30) %>%
  mutate(C_short = case_when(
    grepl("TCGA", Cohort) ~ "TCGA",
    grepl("ICGC", Cohort) ~ "ICGC",
    TRUE                  ~ "GSE14520"
  )) %>%
  select(Signature, C_short, Std_Drift) %>%
  pivot_wider(names_from = C_short, values_from = Std_Drift)

z_tcga_icgc_s <- cor(wide30_zs$TCGA, wide30_zs$ICGC, method = "spearman")
z_tcga_gse_s  <- cor(wide30_zs$TCGA, wide30_zs$GSE14520, method = "spearman")

message("\n>>> [5/6] 跨队列 30% 混入敏感性排序一致性 (Spearman rho):")
message(sprintf("    [ssGSEA] TCGA vs ICGC: rho = %.3f | TCGA vs GSE14520: rho = %.3f | ICGC vs GSE14520: rho = %.3f",
                r_tcga_icgc_s, r_tcga_gse_s, r_icgc_gse_s))
message(sprintf("    [Zscore] TCGA vs ICGC: rho = %.3f | TCGA vs GSE14520: rho = %.3f",
                z_tcga_icgc_s, z_tcga_gse_s))

message("\n>>> 三队列 30% 混入 (ssGSEA) 漂移与翻转率对照表:")
print(as.data.frame(
  wide30_ss %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(Signature,
           TCGA_Drift = Std_Drift_TCGA, ICGC_Drift = Std_Drift_ICGC, GSE_Drift = Std_Drift_GSE14520,
           TCGA_Flip = Flip_Fixed_Pct_TCGA, ICGC_Flip = Flip_Fixed_Pct_ICGC, GSE_Flip = Flip_Fixed_Pct_GSE14520)
), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. 绘制三队列一致性总图 (Fig2C_MultiCohort_Consistency)
# ------------------------------------------------------------------------------
sig_ord <- rev(wide30_ss$Signature)

p_heat <- sim_all %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  mutate(
    Signature = factor(Signature, levels = sig_ord),
    Cohort_lbl = factor(Cohort, levels = unique(sim_all$Cohort))
  ) %>%
  ggplot(aes(x = Cohort_lbl, y = Signature, fill = Std_Drift)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Std_Drift)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Std Drift (30%)") +
  labs(title = "A. Cross-Cohort Drift Landscape (ssGSEA, 30%)", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, face = "bold"),
        plot.title = element_text(face = "bold"))

p_sc1 <- ggplot(wide30_ss, aes(x = Std_Drift_TCGA, y = Std_Drift_ICGC, color = Axis)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7, fill = "grey85") +
  geom_point(size = 2.8) +
  annotate("text", x = min(wide30_ss$Std_Drift_TCGA), y = max(wide30_ss$Std_Drift_ICGC),
           hjust = 0, vjust = 1, fontface = "bold", size = 3.6,
           label = sprintf("Spearman rho = %.3f\n(RNA-seq vs RNA-seq)", r_tcga_icgc_s)) +
  labs(title = "B. TCGA-LIHC vs. ICGC-LIRI-JP",
       x = "TCGA-LIHC Drift (30%)", y = "ICGC-LIRI-JP Drift (30%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

p_sc2 <- ggplot(wide30_ss, aes(x = Std_Drift_TCGA, y = Std_Drift_GSE14520, color = Axis)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7, fill = "grey85") +
  geom_point(size = 2.8) +
  annotate("text", x = min(wide30_ss$Std_Drift_TCGA), y = max(wide30_ss$Std_Drift_GSE14520),
           hjust = 0, vjust = 1, fontface = "bold", size = 3.6,
           label = sprintf("Spearman rho = %.3f\n(RNA-seq vs Microarray)", r_tcga_gse_s)) +
  labs(title = "C. TCGA-LIHC vs. GSE14520",
       x = "TCGA-LIHC Drift (30%)", y = "GSE14520 Drift (30%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p_comb <- (p_heat | (p_sc1 / p_sc2)) + plot_layout(widths = c(1.15, 1))

ggsave("04_figures/Fig2C_MultiCohort_Consistency.pdf", p_comb, width = 13, height = 8)
ggsave("04_figures/Fig2C_MultiCohort_Consistency.png", p_comb, width = 13, height = 8, dpi = 300)

message("\n>>> [6/6] 全部完成！请查看 04_figures/Fig2C_MultiCohort_Consistency.png。")
