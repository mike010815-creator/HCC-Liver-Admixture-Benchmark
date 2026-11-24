# ==============================================================================
# 脚本名称: 02c_gtex_vs_adjacent.R
# 核心任务:
#   1. 补全免疫特征库 (加入 T_cell_exhaustion / Treg / G2M 增殖阴性对照)
#   2. 从 Toil 联合队列中按列轻量提取 TCGA-LIHC 配对样本 + GTEx 健康肝 (零批次差异)
#   3. 头对头对比: 混入“TCGA 癌旁肝” vs 混入“GTEx 健康正常肝” (ssGSEA + Zscore)
#   4. 探测 HCCDB18 与 GSE14520 表头格式，为下一步外部验证铺路
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GSVA)
  library(msigdbr)
  library(IOBR)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 补全免疫特征库并更新 immune_signatures.rds (21 -> 24 个特征)
# ------------------------------------------------------------------------------
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
sig_list <- sigs_obj$sig_list
sig_meta <- sigs_obj$sig_meta

iobr_pool <- c(IOBR::load_data("signature_tme"), IOBR::load_data("signature_collection"))
msig_h <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)

sig_list$T_cell_exhaustion   <- unique(iobr_pool[["T_cell_exhaustion_Peng_et_al"]])
sig_list$Treg_Signature      <- unique(iobr_pool[["Treg_Bindea_et_al"]])
sig_list$CTRL_G2M_Checkpoint <- unique(msig_h$gene_symbol[msig_h$gs_name == "HALLMARK_G2M_CHECKPOINT"])

extra_meta <- tibble(
  Signature     = c("T_cell_exhaustion", "Treg_Signature", "CTRL_G2M_Checkpoint"),
  Axis          = c("2_Suppression_Myeloid", "2_Suppression_Myeloid", "0_Tumor_Proliferation_Ctrl"),
  Total_Genes   = c(length(sig_list$T_cell_exhaustion),
                    length(sig_list$Treg_Signature),
                    length(sig_list$CTRL_G2M_Checkpoint)),
  Detected_TCGA = NA_integer_,
  Detect_Rate   = NA_real_
)

sig_meta <- bind_rows(
  sig_meta %>% filter(!Signature %in% extra_meta$Signature),
  extra_meta
)

saveRDS(list(sig_list = sig_list, sig_meta = sig_meta),
        file = "01_clean_data/immune_signatures.rds")
message(sprintf(">>> [1/5] 免疫特征库已更新，当前共 %d 个特征（含 1 个增殖阴性对照）。",
                length(sig_list)))

# ------------------------------------------------------------------------------
# 2. 解析 Toil 表型表，精准锁定 TCGA-LIHC 配对样本与 GTEx 健康肝样本
# ------------------------------------------------------------------------------
f_toil_tpm   <- "00_raw_data/GTEx_Liver/TcgaTargetGtex_rsem_gene_tpm.gz"
f_toil_pheno <- "00_raw_data/GTEx_Liver/TcgaTargetGTEX_phenotype.txt.gz"
f_toil_probe <- "00_raw_data/GTEx_Liver/gencode.v23.annotation.gene.probemap"

stopifnot(file.exists(f_toil_tpm), file.exists(f_toil_pheno), file.exists(f_toil_probe))

pheno <- fread(f_toil_pheno, data.table = FALSE)
colnames(pheno)[1] <- "sample_id"

col_study   <- grep("study",   colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
col_site    <- grep("primary_site", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
col_disease <- grep("primary disease", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
col_type    <- grep("sample_type", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]

# 仅读取大矩阵第 1 行表头，确认哪些样本真实存在于表达矩阵中
header_cols <- names(fread(f_toil_tpm, nrows = 0))
id_col_name <- header_cols[1]

pheno_in <- pheno %>% filter(sample_id %in% header_cols)

# 筛选 GTEx 健康正常肝
gtex_liver <- pheno_in %>%
  filter(.data[[col_study]] == "GTEX",
         grepl("^Liver$", .data[[col_site]], ignore.case = TRUE))

# 筛选 TCGA 肝细胞癌 (严格排除胆管癌 CHOL)
tcga_lihc <- pheno_in %>%
  filter(.data[[col_study]] == "TCGA",
         grepl("Hepatocellular", .data[[col_disease]], ignore.case = TRUE)) %>%
  mutate(
    patient_id = substr(sample_id, 1, 12),
    type_code  = substr(sample_id, 14, 15)
  )

lihc_T_df <- tcga_lihc %>% filter(type_code == "01") %>% distinct(patient_id, .keep_all = TRUE)
lihc_N_df <- tcga_lihc %>% filter(type_code == "11") %>% distinct(patient_id, .keep_all = TRUE)

paired_pts <- sort(intersect(lihc_T_df$patient_id, lihc_N_df$patient_id))
paired_T   <- lihc_T_df %>% filter(patient_id %in% paired_pts) %>% arrange(patient_id)
paired_N   <- lihc_N_df %>% filter(patient_id %in% paired_pts) %>% arrange(patient_id)

stopifnot(identical(paired_T$patient_id, paired_N$patient_id))
message(sprintf(">>> [2/5] Toil 队列匹配成功: TCGA-LIHC 配对 = %d 对 | GTEx 健康肝 = %d 例",
                nrow(paired_T), nrow(gtex_liver)))

# ------------------------------------------------------------------------------
# 3. 按列轻量读取 1.32 GB 大文件（只读需要的列），按 2^x - 0.001 还原线性 TPM
# ------------------------------------------------------------------------------
cols_need <- unique(c(id_col_name, paired_T$sample_id, paired_N$sample_id, gtex_liver$sample_id))
message(sprintf(">>> [3/5] 正在按列提取 %d 列（跳过其余 18,000+ 列），请稍候约 20 秒...",
                length(cols_need)))

toil_sub <- fread(f_toil_tpm, select = cols_need, data.table = FALSE)
probemap <- fread(f_toil_probe, data.table = FALSE)

ens_ids <- toil_sub[[1]]
mat_log <- as.matrix(toil_sub[, -1, drop = FALSE])
mode(mat_log) <- "numeric"
rownames(mat_log) <- ens_ids

# Toil RSEM TPM 的官方存储单位是 log2(TPM + 0.001)
mat_lin <- 2^mat_log - 0.001
mat_lin[mat_lin < 0] <- 0

# Ensembl (v23) -> Gene Symbol 去重
idx <- match(rownames(mat_lin), probemap$id)
if (sum(!is.na(idx)) < 1000) {
  idx <- match(sub("\\..*$", "", rownames(mat_lin)), sub("\\..*$", "", probemap$id))
}
gsym <- probemap$gene[idx]
keep <- !is.na(gsym) & gsym != "" & !grepl("^\\s*$", gsym)
mat_lin <- mat_lin[keep, , drop = FALSE]
gsym    <- gsym[keep]

ord <- order(gsym, -rowMeans(mat_lin, na.rm = TRUE))
mat_lin <- mat_lin[ord, , drop = FALSE]
gsym    <- gsym[ord]
mat_dedup <- mat_lin[!duplicated(gsym), , drop = FALSE]
rownames(mat_dedup) <- gsym[!duplicated(gsym)]
mat_dedup <- mat_dedup[rowSums(mat_dedup) > 0, , drop = FALSE]

# 每列重标化至严格 1e6
mat_norm <- sweep(mat_dedup, 2, colSums(mat_dedup), FUN = "/") * 1e6

T_toil    <- mat_norm[, paired_T$sample_id,   drop = FALSE]
N_adj     <- mat_norm[, paired_N$sample_id,   drop = FALSE]
N_gtex    <- mat_norm[, gtex_liver$sample_id, drop = FALSE]

saveRDS(
  list(T_toil = T_toil, N_adj = N_adj, N_gtex = N_gtex,
       paired_pts = paired_pts),
  file = "01_clean_data/toil_lihc_gtex_clean.rds"
)
message(sprintf("    已保存轻量化子集至 01_clean_data/toil_lihc_gtex_clean.rds (%d 基因)",
                nrow(mat_norm)))

# ------------------------------------------------------------------------------
# 4. 头对头模拟: 30% TCGA 癌旁混入 vs 30% GTEx 健康肝混入 (5 次随机配对取均值)
# ------------------------------------------------------------------------------
run_ssgsea <- function(mat_log2, gsets) {
  p <- ssgseaParam(as.matrix(mat_log2), gsets, minSize = 2, normalize = FALSE)
  gsva(p, verbose = FALSE)
}

base_log <- log2(T_toil + 1)
base_mu  <- rowMeans(base_log)
base_sd  <- apply(base_log, 1, sd)
base_sd[base_sd == 0] <- NA

run_zscore <- function(mat_log2, gsets) {
  zm <- sweep(sweep(mat_log2, 1, base_mu, "-"), 1, base_sd, "/")
  out <- sapply(gsets, function(g) {
    gu <- intersect(g, rownames(zm))
    colMeans(zm[gu, , drop = FALSE], na.rm = TRUE)
  })
  t(out)
}

eval_one_mix <- function(N_mat, alpha = 0.30) {
  mix_log <- log2((1 - alpha) * T_toil + alpha * N_mat + 1)
  s_ss <- run_ssgsea(mix_log, sig_list)
  s_zs <- run_zscore(mix_log, sig_list)
  
  calc_one <- function(s0, sa, method_lbl) {
    map_dfr(rownames(s0), function(sig) {
      v0 <- as.numeric(s0[sig, ]); va <- as.numeric(sa[sig, ])
      med0 <- median(v0)
      tibble(
        Method    = method_lbl,
        Signature = sig,
        Std_Drift = mean(va - v0) / sd(v0),
        Rho       = cor(v0, va, method = "spearman"),
        Flip_Pct  = mean((v0 >= med0) != (va >= med0)) * 100
      )
    })
  }
  bind_rows(calc_one(s0_ss, s_ss, "ssGSEA"), calc_one(s0_zs, s_zs, "Zscore"))
}

message(">>> [4/5] 正在运行 30% 癌旁肝 vs 30% GTEx 健康肝头对头对比...")
s0_ss <- run_ssgsea(base_log, sig_list)
s0_zs <- run_zscore(base_log, sig_list)

# 癌旁混入
res_adj <- eval_one_mix(N_adj, alpha = 0.30) %>% mutate(Source = "TCGA_Adjacent_Liver")

# GTEx 健康肝随机抽样配对 5 次取均值，消除单次抽样偶然性
n_pairs <- ncol(T_toil)
res_gtex <- map_dfr(1:5, function(rep_i) {
  set.seed(2026 + rep_i)
  pick <- sample(ncol(N_gtex), n_pairs, replace = FALSE)
  eval_one_mix(N_gtex[, pick, drop = FALSE], alpha = 0.30)
}) %>%
  group_by(Method, Signature) %>%
  summarise(
    Std_Drift = mean(Std_Drift),
    Rho       = mean(Rho),
    Flip_Pct  = mean(Flip_Pct),
    .groups   = "drop"
  ) %>%
  mutate(Source = "GTEx_Healthy_Liver")

cmp_df <- bind_rows(res_adj, res_gtex) %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature")

cmp_wide <- cmp_df %>%
  select(Method, Signature, Axis, Source, Std_Drift, Flip_Pct) %>%
  pivot_wider(names_from = Source, values_from = c(Std_Drift, Flip_Pct)) %>%
  mutate(Delta_Adj_minus_GTEx = Std_Drift_TCGA_Adjacent_Liver - Std_Drift_GTEx_Healthy_Liver)

message("\n>>> [结果预览 A] Zscore (固定基线) 下 30% 癌旁肝 vs 30% GTEx 健康肝漂移对比:")
print(as.data.frame(
  cmp_wide %>%
    filter(Method == "Zscore") %>%
    arrange(desc(Delta_Adj_minus_GTEx)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(Signature,
           Adj_Drift = Std_Drift_TCGA_Adjacent_Liver,
           GTEx_Drift = Std_Drift_GTEx_Healthy_Liver,
           Diff = Delta_Adj_minus_GTEx,
           Adj_Flip = Flip_Pct_TCGA_Adjacent_Liver,
           GTEx_Flip = Flip_Pct_GTEx_Healthy_Liver)
), row.names = FALSE)

message("\n>>> [结果预览 B] ssGSEA 下 30% 癌旁肝 vs 30% GTEx 健康肝漂移对比:")
print(as.data.frame(
  cmp_wide %>%
    filter(Method == "ssGSEA") %>%
    arrange(desc(Delta_Adj_minus_GTEx)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(Signature,
           Adj_Drift = Std_Drift_TCGA_Adjacent_Liver,
           GTEx_Drift = Std_Drift_GTEx_Healthy_Liver,
           Diff = Delta_Adj_minus_GTEx,
           Adj_Flip = Flip_Pct_TCGA_Adjacent_Liver,
           GTEx_Flip = Flip_Pct_GTEx_Healthy_Liver)
), row.names = FALSE)

write.csv(cmp_df, "03_results/02c_adjacent_vs_gtex_30pct.csv", row.names = FALSE)

# 哑铃对比图 (同时展示 Zscore 和 ssGSEA)
sig_ord <- cmp_wide %>%
  filter(Method == "Zscore") %>%
  arrange(Std_Drift_TCGA_Adjacent_Liver) %>%
  pull(Signature)

p_cmp <- cmp_df %>%
  mutate(Signature = factor(Signature, levels = sig_ord)) %>%
  ggplot(aes(x = Std_Drift, y = Signature)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(aes(group = Signature), color = "grey65", linewidth = 0.9) +
  geom_point(aes(color = Source), size = 3) +
  facet_wrap(~ Method, scales = "free_x") +
  scale_color_manual(
    values = c("TCGA_Adjacent_Liver" = "#D62728", "GTEx_Healthy_Liver" = "#1F77B4"),
    labels = c("GTEx Healthy Liver (30%)", "TCGA Adjacent Liver (30%)"),
    name   = "Admixture Source"
  ) +
  labs(
    title    = "Adjacent Cirrhotic/Inflamed Liver vs. Healthy GTEx Liver Admixture",
    subtitle = "Toil-harmonized cohort (paired TCGA-LIHC tumors, 30% linear TPM admixture)",
    x        = "Standardized Score Drift (in baseline SD units)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F2F4F7"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave("04_figures/Fig2B_Adjacent_vs_GTEx_Healthy.pdf", p_cmp, width = 11, height = 7.5)
ggsave("04_figures/Fig2B_Adjacent_vs_GTEx_Healthy.png", p_cmp, width = 11, height = 7.5, dpi = 300)

# ------------------------------------------------------------------------------
# 5. 顺手探测 HCCDB18 与 GSE14520 表头（确保下一步外部验证脚本零报错）
# ------------------------------------------------------------------------------
message("\n>>> [5/5] 外部验证队列原始文件表头预览:")
cat("\n--- HCCDB18.sample.txt (前 3 行) ---\n")
print(head(fread("00_raw_data/HCCDB18_ICGC/HCCDB18.sample.txt", data.table = FALSE), 3))

cat("\n--- HCCDB18.patient.txt (前 3 行) ---\n")
print(head(fread("00_raw_data/HCCDB18_ICGC/HCCDB18.patient.txt", data.table = FALSE), 3))

cat("\n--- HCCDB18_mRNA_level3.txt (前 3 行 x 前 5 列) ---\n")
h18_peek <- fread("00_raw_data/HCCDB18_ICGC/HCCDB18_mRNA_level3.txt", nrows = 3, data.table = FALSE)
print(h18_peek[, 1:min(5, ncol(h18_peek))])

cat("\n--- GSE14520_Extra_Supplement.txt.gz (前 3 行 x 前 10 列) ---\n")
gse_sup <- fread("00_raw_data/GSE14520/GSE14520_Extra_Supplement.txt.gz", nrows = 3, data.table = FALSE)
print(gse_sup[, 1:min(10, ncol(gse_sup))])

message("\n>>> 全部完成！请查看 04_figures/Fig2B_Adjacent_vs_GTEx_Healthy.png。")
