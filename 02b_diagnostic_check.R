# ==============================================================================
# 脚本名称: 02b_diagnostic_check.R
# 目的: 5 秒排雷自检 —— 区分“真实生物学抬高”与“ssGSEA 排序数学效应”
# ==============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(msigdbr)
  library(IOBR)
})

tcga <- readRDS("01_clean_data/tcga_lihc_clean.rds")
sigs <- readRDS("01_clean_data/immune_signatures.rds")
res2 <- readRDS("01_clean_data/02_tcga_admixture_results.rds")

# ------------------------------------------------------------------------------
# 检查 1: 对比 30% 混入时 ssGSEA / GSVA / Zscore 三种算法的漂移方向 (Std_Drift)
# ------------------------------------------------------------------------------
drift_cmp <- res2$metrics_long %>%
  filter(Alpha_Pct == 30) %>%
  select(Signature, Method, Std_Drift) %>%
  pivot_wider(names_from = Method, values_from = Std_Drift) %>%
  arrange(desc(Zscore))

message(">>> [检查 1] 30% 混入时三种算法的标准化漂移方向对比 (重点看 Zscore 有没有负数):")
print(as.data.frame(drift_cmp %>% mutate(across(where(is.numeric), ~ round(.x, 2)))),
      row.names = FALSE)

# ------------------------------------------------------------------------------
# 检查 2: 加入 3 个肿瘤增殖对照特征，验证本流程能正常检出“负向稀释塌陷”
# ------------------------------------------------------------------------------
msig_h <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
ctrl_sigs <- list(
  CTRL_G2M_Checkpoint = unique(msig_h$gene_symbol[msig_h$gs_name == "HALLMARK_G2M_CHECKPOINT"]),
  CTRL_E2F_Targets    = unique(msig_h$gene_symbol[msig_h$gs_name == "HALLMARK_E2F_TARGETS"]),
  CTRL_MYC_Targets_V1 = unique(msig_h$gene_symbol[msig_h$gs_name == "HALLMARK_MYC_TARGETS_V1"])
)

m0  <- log2(tcga$T_tpm + 1)
m30 <- log2(0.7 * tcga$T_tpm + 0.3 * tcga$N_tpm + 1)
m100 <- log2(tcga$N_tpm + 1)

run_ssgsea_quick <- function(mat, gsets) {
  p <- ssgseaParam(as.matrix(mat), gsets, minSize = 2, normalize = FALSE)
  gsva(p, verbose = FALSE)
}
s0   <- run_ssgsea_quick(m0, ctrl_sigs)
s30  <- run_ssgsea_quick(m30, ctrl_sigs)
s100 <- run_ssgsea_quick(m100, ctrl_sigs)

# 同时看 Zscore 漂移
mu0 <- rowMeans(m0); sd0 <- apply(m0, 1, sd); sd0[sd0 == 0] <- NA
z_drift <- sapply(ctrl_sigs, function(g) {
  gu <- intersect(g, rownames(m0))
  z0  <- colMeans(sweep(sweep(m0[gu, ],  1, mu0[gu], "-"), 1, sd0[gu], "/"), na.rm = TRUE)
  z30 <- colMeans(sweep(sweep(m30[gu, ], 1, mu0[gu], "-"), 1, sd0[gu], "/"), na.rm = TRUE)
  mean(z30 - z0) / sd(z0)
})

ctrl_tbl <- tibble(
  Control_Signature  = names(ctrl_sigs),
  ssGSEA_Drift_30pct = round(sapply(names(ctrl_sigs), function(k) mean(s30[k, ] - s0[k, ]) / sd(s0[k, ])), 2),
  ssGSEA_Drift_100pct = round(sapply(names(ctrl_sigs), function(k) mean(s100[k, ] - s0[k, ]) / sd(s0[k, ])), 2),
  Zscore_Drift_30pct = round(as.numeric(z_drift), 2)
)
message("\n>>> [检查 2] 肿瘤增殖对照特征的漂移方向 (预期应为明显负数):")
print(as.data.frame(ctrl_tbl), row.names = FALSE)

# ------------------------------------------------------------------------------
# 检查 3: 拆解 WNT_target 里到底包含哪些基因，为何飙升到 +0.83 SD
# ------------------------------------------------------------------------------
wnt_genes <- intersect(sigs$sig_list$WNT_target, rownames(tcga$T_tpm))
wnt_tbl <- tibble(
  Gene = wnt_genes,
  Tumor_Med_TPM  = round(apply(tcga$T_tpm[wnt_genes, , drop = FALSE], 1, median), 1),
  Normal_Med_TPM = round(apply(tcga$N_tpm[wnt_genes, , drop = FALSE], 1, median), 1),
  Ratio_N_over_T = round((Normal_Med_TPM + 0.1) / (Tumor_Med_TPM + 0.1), 2)
) %>% arrange(desc(Ratio_N_over_T))
message("\n>>> [检查 3] WNT_target 组成基因在癌 vs 癌旁中的中位 TPM:")
print(as.data.frame(wnt_tbl), row.names = FALSE)

# ------------------------------------------------------------------------------
# 检查 4: 顺手查本机 IOBR 中 Exhaustion 与 Treg 的准确键名
# ------------------------------------------------------------------------------
pool <- c(IOBR::load_data("signature_tme"), IOBR::load_data("signature_collection"))
message("\n>>> [检查 4] 本机 IOBR 中含 exhaust / treg 的可用特征名:")
print(grep("exhaust|treg", names(pool), ignore.case = TRUE, value = TRUE))
