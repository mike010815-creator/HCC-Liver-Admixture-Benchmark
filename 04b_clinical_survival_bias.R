# ==============================================================================
# 脚本名称: 04b_clinical_survival_bias.R
# 核心任务:
#   1. 整合三大初治 HCC 临床队列 (TCGA-LIHC, ICGC-LIRI-JP, GSE14520, ~790 例原发癌)
#   2. 真实世界共变验证: 检验真实未混合肿瘤中 Spearman(NHAS, Signature) 是否复现
#      in silico 30% 模拟敏感性排序 (且增殖对照 CTRL_G2M 呈负相关)
#   3. 真实预后偏倚验证: Low-NHAS vs High-NHAS 分层 Cox 回归 + NHAS 校正前后 HR 漂移
#      + 三队列 Meta 汇总，导出 Fig5_RealWorld_Clinical_Bias
#   4. 3 秒顺手诊断 GSE238264 (HCC3R) 的 UMI 混杂与原作者病理注释
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Matrix)
  library(GSVA)
  library(survival)
  library(Seurat)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 加载三大清洗好的 Bulk 队列、免疫特征库、NHAS 基因与模拟基准结果
# ------------------------------------------------------------------------------
tcga_obj <- readRDS("01_clean_data/tcga_lihc_clean.rds")
h18_obj  <- readRDS("01_clean_data/hccdb18_clean.rds")
gse_obj  <- readRDS("01_clean_data/gse14520_clean.rds")
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
bulk_sim <- readRDS("01_clean_data/02d_multicohort_sim_results.rds")

sig_list   <- sigs_obj$sig_list
sig_meta   <- sigs_obj$sig_meta
nhas_genes <- res3_obj$nhas_genes

all_gsets <- c(list(NHAS = nhas_genes), sig_list)

# ------------------------------------------------------------------------------
# 2. 统一整理三大队列原发癌的 OS 生存表 (统一时间单位为月, OS_event 为 0/1)
# ------------------------------------------------------------------------------
# 2.1 TCGA-LIHC (OS.time 为天 -> 转月)
tcga_surv <- tcga_obj$all_T_meta %>%
  transmute(
    sample_id = sample_id,
    OS_time   = suppressWarnings(as.numeric(OS.time)) / 30.4375,
    OS_event  = suppressWarnings(as.numeric(OS))
  ) %>%
  filter(sample_id %in% colnames(tcga_obj$all_T_tpm),
         !is.na(OS_time), !is.na(OS_event), OS_time > 0)

# 2.2 ICGC-LIRI-JP (HCCDB18): 自动匹配样本表与患者表并解析生存字段
h18_t_samp <- h18_obj$samp_meta %>%
  filter(sample_id %in% colnames(h18_obj$all_T_tpm))
h18_pat <- h18_obj$pat_meta

# 自动选择交集最大的患者关联键
ov_pat_id <- sum(h18_t_samp$pat_key %in% h18_pat$patient_id)
ov_pat_do <- if ("PATIENT" %in% colnames(h18_t_samp) && "PATIENT" %in% colnames(h18_pat)) {
  sum(h18_t_samp$PATIENT %in% h18_pat$PATIENT)
} else 0

if (ov_pat_id >= ov_pat_do) {
  h18_merged <- h18_t_samp %>% left_join(h18_pat, by = c("pat_key" = "patient_id"))
} else {
  h18_merged <- h18_t_samp %>% left_join(h18_pat, by = "PATIENT")
}

# 解析 STATUS (兼容 Alive/Dead 或 0/1) 与 SUR (自动判断天/月)
raw_stat <- trimws(as.character(h18_merged$STATUS))
raw_sur  <- suppressWarnings(as.numeric(h18_merged$SUR))
message(">>> [1/5] HCCDB18 原始生存字段核对:")
message("    STATUS 取值分布: ", paste(names(table(raw_stat)), table(raw_stat), sep = "=", collapse = " | "))
message(sprintf("    SUR 数值中位数: %.1f | 最大值: %.1f", median(raw_sur, na.rm = TRUE), max(raw_sur, na.rm = TRUE)))

h18_event <- case_when(
  grepl("Dead|Deceased|^1$", raw_stat, ignore.case = TRUE) ~ 1,
  grepl("Alive|Living|^0$",  raw_stat, ignore.case = TRUE) ~ 0,
  TRUE ~ NA_real_
)
# 若 SUR 最大值 > 300 视为天数转月，否则视为月数
h18_time_mo <- if (max(raw_sur, na.rm = TRUE) > 300) raw_sur / 30.4375 else raw_sur

h18_surv <- tibble(
  sample_id = h18_merged$sample_id,
  OS_time   = h18_time_mo,
  OS_event  = h18_event
) %>%
  filter(!is.na(OS_time), !is.na(OS_event), OS_time > 0)

# 2.3 GSE14520 (GPL3921): Survival months & Survival status
gse_meta <- gse_obj$meta %>% filter(sample_id %in% colnames(gse_obj$all_T_lin))
col_gse_t <- grep("Survival.*month", colnames(gse_meta), ignore.case = TRUE, value = TRUE)[1]
col_gse_e <- grep("Survival.*status", colnames(gse_meta), ignore.case = TRUE, value = TRUE)[1]

gse_surv <- tibble(
  sample_id = gse_meta$sample_id,
  OS_time   = suppressWarnings(as.numeric(gse_meta[[col_gse_t]])),
  OS_event  = suppressWarnings(as.numeric(gse_meta[[col_gse_e]]))
) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  filter(!is.na(OS_time), !is.na(OS_event), OS_time > 0)

message(sprintf("    有效总生存队列汇总: TCGA-LIHC (n=%d, 死亡=%d) | ICGC-LIRI-JP (n=%d, 死亡=%d) | GSE14520 (n=%d, 死亡=%d)",
                nrow(tcga_surv), sum(tcga_surv$OS_event == 1),
                nrow(h18_surv),  sum(h18_surv$OS_event == 1),
                nrow(gse_surv),  sum(gse_surv$OS_event == 1)))

# ------------------------------------------------------------------------------
# 3. 对三大真实队列全部原发癌跑单样本 ssGSEA (NHAS + 24 个特征)
# ------------------------------------------------------------------------------
score_cohort_real <- function(lin_mat, surv_df, cohort_lbl) {
  common_s <- intersect(colnames(lin_mat), surv_df$sample_id)
  m_log <- log2(lin_mat[, common_s, drop = FALSE] + 1)
  p <- ssgseaParam(as.matrix(m_log), all_gsets, minSize = 2, normalize = FALSE)
  sc_mat <- gsva(p, verbose = FALSE)
  
  sc_df <- as.data.frame(t(sc_mat)) %>%
    rownames_to_column("sample_id") %>%
    inner_join(surv_df, by = "sample_id") %>%
    mutate(Cohort = cohort_lbl)
  return(sc_df)
}

message("\n>>> [2/5] 正在三大真实原发肝癌队列上计算 NHAS 与 24 个特征得分...")
real_tcga <- score_cohort_real(tcga_obj$all_T_tpm, tcga_surv, "TCGA_LIHC")
real_icgc <- score_cohort_real(h18_obj$all_T_tpm,  h18_surv,  "ICGC_LIRI_JP")
real_gse  <- score_cohort_real(gse_obj$all_T_lin,  gse_surv,  "GSE14520")

real_all_list <- list(TCGA_LIHC = real_tcga, ICGC_LIRI_JP = real_icgc, GSE14520 = real_gse)
sig_names <- names(sig_list)

# ------------------------------------------------------------------------------
# 4. 真实世界验证 1: 未混合真实原发癌中 Spearman(NHAS, Signature) vs Bulk 模拟漂移
# ------------------------------------------------------------------------------
real_cor_df <- map_dfr(names(real_all_list), function(cn) {
  df <- real_all_list[[cn]]
  map_dfr(sig_names, function(s) {
    ct <- suppressWarnings(cor.test(df$NHAS, df[[s]], method = "spearman"))
    tibble(
      Cohort    = cn,
      Signature = s,
      Real_Rho_with_NHAS = as.numeric(ct$estimate),
      Real_P             = as.numeric(ct$p.value)
    )
  })
})

real_cor_wide <- real_cor_df %>%
  select(Cohort, Signature, Real_Rho_with_NHAS) %>%
  pivot_wider(names_from = Cohort, values_from = Real_Rho_with_NHAS, names_prefix = "RealRho_") %>%
  mutate(Mean_Real_Rho = (RealRho_TCGA_LIHC + RealRho_ICGC_LIRI_JP + RealRho_GSE14520) / 3) %>%
  left_join(
    bulk_sim %>%
      filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
      mutate(C_short = case_when(
        grepl("TCGA", Cohort) ~ "TCGA",
        grepl("ICGC", Cohort) ~ "ICGC",
        TRUE                  ~ "GSE14520"
      )) %>%
      select(Signature, C_short, Std_Drift) %>%
      pivot_wider(names_from = C_short, values_from = Std_Drift, names_prefix = "Sim30_"),
    by = "Signature"
  ) %>%
  mutate(Mean_Sim30_Drift = (Sim30_TCGA + Sim30_ICGC + Sim30_GSE14520) / 3) %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature") %>%
  arrange(desc(Mean_Real_Rho))

c_tcga <- cor(real_cor_wide$RealRho_TCGA_LIHC,    real_cor_wide$Sim30_TCGA,     method = "spearman")
c_icgc <- cor(real_cor_wide$RealRho_ICGC_LIRI_JP, real_cor_wide$Sim30_ICGC,     method = "spearman")
c_gse  <- cor(real_cor_wide$RealRho_GSE14520,     real_cor_wide$Sim30_GSE14520, method = "spearman")
c_mean <- cor(real_cor_wide$Mean_Real_Rho,        real_cor_wide$Mean_Sim30_Drift, method = "spearman")

message("\n>>> [3/5] 真实原发癌中 Spearman(NHAS, Signature) vs 30% 模拟敏感性一致性:")
message(sprintf("    TCGA-LIHC (n=%d): rho = %.3f | ICGC-LIRI-JP (n=%d): rho = %.3f | GSE14520 (n=%d): rho = %.3f | 三队列均值: rho = %.3f",
                nrow(real_tcga), c_tcga, nrow(real_icgc), c_icgc, nrow(real_gse), c_gse, c_mean))

print(as.data.frame(
  real_cor_wide %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    select(Signature, RealRho_TCGA = RealRho_TCGA_LIHC, RealRho_ICGC = RealRho_ICGC_LIRI_JP,
           RealRho_GSE = RealRho_GSE14520, Mean_Real_Rho, Mean_Sim30_Drift)
), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. 真实世界验证 2: 生存预后偏倚分析 (Low-NHAS vs High-NHAS 分层 Cox + NHAS 校正 + Meta)
# ------------------------------------------------------------------------------
fit_cox_one <- function(df, formula_str, term_name) {
  fit <- tryCatch(coxph(as.formula(formula_str), data = df), error = function(e) NULL)
  if (is.null(fit)) return(c(beta = NA, se = NA, p = NA))
  sm <- summary(fit)$coefficients
  if (!(term_name %in% rownames(sm))) return(c(beta = NA, se = NA, p = NA))
  c(beta = sm[term_name, "coef"], se = sm[term_name, "se(coef)"], p = sm[term_name, "Pr(>|z|)"])
}

cox_cohort_df <- map_dfr(names(real_all_list), function(cn) {
  df <- real_all_list[[cn]]
  df$NHAS_z <- as.numeric(scale(df$NHAS))
  med_nhas  <- median(df$NHAS)
  df_low    <- df %>% filter(NHAS <  med_nhas)
  df_high   <- df %>% filter(NHAS >= med_nhas)
  
  # 先看 NHAS 本身对 OS 的预后效应 (供审稿背景说明)
  nhas_os <- fit_cox_one(df, "Surv(OS_time, OS_event) ~ NHAS_z", "NHAS_z")
  message(sprintf("    [%s] NHAS 本身与总生存 OS 关联: HR = %.2f (P = %.3g)",
                  cn, exp(nhas_os["beta"]), nhas_os["p"]))
  
  map_dfr(sig_names, function(s) {
    df$Sig_z      <- as.numeric(scale(df[[s]]))
    df_low$Sig_z  <- as.numeric(scale(df_low[[s]]))
    df_high$Sig_z <- as.numeric(scale(df_high[[s]]))
    
    c_raw  <- fit_cox_one(df,      "Surv(OS_time, OS_event) ~ Sig_z", "Sig_z")
    c_adj  <- fit_cox_one(df,      "Surv(OS_time, OS_event) ~ Sig_z + NHAS_z", "Sig_z")
    c_low  <- fit_cox_one(df_low,  "Surv(OS_time, OS_event) ~ Sig_z", "Sig_z")
    c_high <- fit_cox_one(df_high, "Surv(OS_time, OS_event) ~ Sig_z", "Sig_z")
    c_int  <- fit_cox_one(df,      "Surv(OS_time, OS_event) ~ Sig_z * NHAS_z", "Sig_z:NHAS_z")
    
    tibble(
      Cohort     = cn,
      Signature  = s,
      beta_raw   = c_raw["beta"],  se_raw  = c_raw["se"],  p_raw  = c_raw["p"],
      beta_adj   = c_adj["beta"],  se_adj  = c_adj["se"],  p_adj  = c_adj["p"],
      beta_low   = c_low["beta"],  se_low  = c_low["se"],  p_low  = c_low["p"],
      beta_high  = c_high["beta"], se_high = c_high["se"], p_high = c_high["p"],
      beta_int   = c_int["beta"],  se_int  = c_int["se"],  p_int  = c_int["p"]
    )
  })
})

# 三队列固定效应逆方差加权 Meta 汇总 (合并近 790 例患者统计检验力)
pool_ivw <- function(beta, se) {
  ok <- !is.na(beta) & !is.na(se) & se > 0
  w  <- 1 / (se[ok]^2)
  b  <- sum(w * beta[ok]) / sum(w)
  s  <- sqrt(1 / sum(w))
  p  <- 2 * pnorm(-abs(b / s))
  list(hr = exp(b), lo = exp(b - 1.96 * s), hi = exp(b + 1.96 * s), beta = b, se = s, p = p)
}

meta_cox_df <- cox_cohort_df %>%
  group_by(Signature) %>%
  summarise(
    raw_res  = list(pool_ivw(beta_raw,  se_raw)),
    adj_res  = list(pool_ivw(beta_adj,  se_adj)),
    low_res  = list(pool_ivw(beta_low,  se_low)),
    high_res = list(pool_ivw(beta_high, se_high)),
    int_res  = list(pool_ivw(beta_int,  se_int)),
    .groups  = "drop"
  ) %>%
  mutate(
    HR_Raw   = sapply(raw_res,  `[[`, "hr"),   P_Raw  = sapply(raw_res,  `[[`, "p"),
    HR_Adj   = sapply(adj_res,  `[[`, "hr"),   P_Adj  = sapply(adj_res,  `[[`, "p"),
    HR_Low   = sapply(low_res,  `[[`, "hr"),   Lo_Low = sapply(low_res,  `[[`, "lo"),  Hi_Low = sapply(low_res,  `[[`, "hi"),  P_Low  = sapply(low_res,  `[[`, "p"),
    HR_High  = sapply(high_res, `[[`, "hr"),   Lo_High = sapply(high_res, `[[`, "lo"), Hi_High = sapply(high_res, `[[`, "hi"), P_High = sapply(high_res, `[[`, "p"),
    Delta_logHR_Adj_minus_Raw = sapply(adj_res, `[[`, "beta") - sapply(raw_res, `[[`, "beta"),
    Delta_logHR_High_minus_Low = sapply(high_res, `[[`, "beta") - sapply(low_res, `[[`, "beta"),
    P_Int_Meta = sapply(int_res, `[[`, "p")
  ) %>%
  select(-ends_with("_res")) %>%
  left_join(real_cor_wide %>% select(Signature, Axis, Mean_Real_Rho, Mean_Sim30_Drift), by = "Signature") %>%
  arrange(desc(abs(Delta_logHR_Adj_minus_Raw)))

message("\n>>> [4/5] 三队列 Meta 汇总 (n~790): NHAS 校正前后及低/高混入分层 HR 变化总表:")
print(as.data.frame(
  meta_cox_df %>%
    mutate(across(c(HR_Raw, HR_Adj, HR_Low, HR_High, Delta_logHR_Adj_minus_Raw, Mean_Real_Rho), ~ round(.x, 3)),
           across(c(P_Raw, P_Adj, P_Low, P_High, P_Int_Meta), ~ signif(.x, 2))) %>%
    select(Signature, Mean_Real_Rho, HR_Raw, P_Raw, HR_Adj, P_Adj, HR_Low, P_Low, HR_High, P_High, P_Int_Meta)
), row.names = FALSE)

# 保存真实临床验证全部结果
saveRDS(
  list(
    real_all_list = real_all_list,
    real_cor_df   = real_cor_df,
    real_cor_wide = real_cor_wide,
    cox_cohort_df = cox_cohort_df,
    meta_cox_df   = meta_cox_df
  ),
  file = "01_clean_data/04b_clinical_survival_results.rds"
)
write.csv(real_cor_wide, "03_results/04b_realworld_nhas_correlation.csv", row.names = FALSE)
write.csv(meta_cox_df,   "03_results/04b_meta_cox_prognostic_bias.csv",   row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. 绘制 Figure 5: 真实肿瘤共变验证 + 生存 HR 混杂偏移总图
# ------------------------------------------------------------------------------
sig_ord5 <- rev(real_cor_wide$Signature)

# 图 5A: 三大真实队列中每个特征与 NHAS 的 Spearman 相关热图
p5a <- real_cor_df %>%
  mutate(
    Signature = factor(Signature, levels = sig_ord5),
    Cohort_lbl = factor(Cohort, levels = c("TCGA_LIHC", "ICGC_LIRI_JP", "GSE14520"),
                        labels = c(sprintf("TCGA-LIHC\n(n=%d)", nrow(real_tcga)),
                                   sprintf("ICGC-LIRI-JP\n(n=%d)", nrow(real_icgc)),
                                   sprintf("GSE14520\n(n=%d)", nrow(real_gse)))),
    sig_star = ifelse(Real_P < 0.001, "***", ifelse(Real_P < 0.01, "**", ifelse(Real_P < 0.05, "*", "")))
  ) %>%
  ggplot(aes(x = Cohort_lbl, y = Signature, fill = Real_Rho_with_NHAS)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f%s", Real_Rho_with_NHAS, sig_star)), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Spearman rho\nwith NHAS") +
  labs(title = "A. Real-World Co-variation with NHAS (Unmixed Tumors)", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(face = "bold"), plot.title = element_text(face = "bold"))

# 图 5B: 真实肿瘤中 Spearman(NHAS, Sig) vs 30% 模拟漂移跨验证散点图
p5b <- ggplot(real_cor_wide, aes(x = Mean_Sim30_Drift, y = Mean_Real_Rho, color = Axis)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_smooth(method = "lm", se = TRUE, color = "grey25", linewidth = 0.8, fill = "grey88") +
  geom_point(size = 3) +
  annotate("text", x = min(real_cor_wide$Mean_Sim30_Drift), y = max(real_cor_wide$Mean_Real_Rho),
           hjust = 0, vjust = 1, fontface = "bold", size = 3.5,
           label = sprintf("Pooled 3-Cohort: rho = %.3f\n(TCGA %.2f | ICGC %.2f | GSE %.2f)",
                           c_mean, c_tcga, c_icgc, c_gse)) +
  labs(
    title    = "B. In Silico Prediction vs. Real-World Observation",
    subtitle = "Mean across 3 treatment-naive HCC cohorts (~790 patients)",
    x        = "Predicted 30% Admixture Drift (Mean Std Drift)",
    y        = "Observed Real-Tumor Correlation with NHAS (Mean rho)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(face = "bold"))

# 图 5C: NHAS 校正前 (Raw) vs 校正后 (Adjusted) 的 Meta 风险比偏移哑铃图
p5c_df <- meta_cox_df %>%
  select(Signature, HR_Raw, HR_Adj, P_Raw, P_Adj) %>%
  pivot_longer(cols = c(HR_Raw, HR_Adj), names_to = "Model", values_to = "HR") %>%
  mutate(
    Signature = factor(Signature, levels = sig_ord5),
    Model = recode(Model, HR_Raw = "Unadjusted Cox (~ Sig)", HR_Adj = "NHAS-Adjusted Cox (~ Sig + NHAS)")
  )

p5c <- ggplot(p5c_df, aes(x = HR, y = Signature)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40") +
  geom_line(aes(group = Signature), color = "grey65", linewidth = 0.9) +
  geom_point(aes(color = Model, shape = Model), size = 3) +
  scale_color_manual(values = c("Unadjusted Cox (~ Sig)" = "#D62728",
                                "NHAS-Adjusted Cox (~ Sig + NHAS)" = "#1F77B4"),
                     name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(
    title    = "C. Prognostic Confounding by NHAS (3-Cohort Meta-Cox)",
    subtitle = "Pooled Hazard Ratio per 1-SD increase in signature score",
    x        = "Pooled Hazard Ratio (OS)", y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "vertical",
        plot.title = element_text(face = "bold"))

p5 <- (p5a | p5b | p5c) + plot_layout(widths = c(1.1, 1.15, 1.05))
ggsave("04_figures/Fig5_RealWorld_Clinical_Bias.pdf", p5, width = 16, height = 8.5)
ggsave("04_figures/Fig5_RealWorld_Clinical_Bias.png", p5, width = 16, height = 8.5, dpi = 300)

# ------------------------------------------------------------------------------
# 7. [附带 3 秒诊断] 检查 GSE238264 (HCC3R) 的 UMI 深度混杂与作者 SME 病理注释
# ------------------------------------------------------------------------------
message("\n>>> [5/5] 附带诊断 GSE238264 (HCC3R) 测序深度混杂与作者聚类注释:")
h3_dir <- "01_clean_data/spatial_unpacked/GSE238264/HCC3R"
if (dir.exists(h3_dir)) {
  cnt3 <- Read10X(file.path(h3_dir, "filtered_feature_bc_matrix"))
  if (is.list(cnt3)) cnt3 <- cnt3[[1]]
  seu3 <- CreateSeuratObject(counts = cnt3, assay = "Spatial")
  seu3 <- subset(seu3, subset = nFeature_Spatial >= 200)
  seu3 <- NormalizeData(seu3, verbose = FALSE)
  
  # 用 AddModuleScore (自动扣除匹配表达量背景对照基因，消除 UMI/dropout 混杂)
  nhas_u <- intersect(nhas_genes, rownames(seu3))
  seu3 <- AddModuleScore(seu3, features = list(nhas_u), name = "NHAS_MS_", ctrl = 50)
  seu3$NHAS_RawMean <- Matrix::colMeans(GetAssayData(seu3, layer = "data")[nhas_u, , drop = FALSE])
  
  g2m_u  <- intersect(sig_list$CTRL_G2M_Checkpoint, rownames(seu3))
  comp_u <- intersect(sig_list$HM_Complement, rownames(seu3))
  tex_u  <- intersect(sig_list$T_cell_exhaustion, rownames(seu3))
  seu3 <- AddModuleScore(seu3, features = list(g2m_u, comp_u, tex_u), name = "SIG_MS_", ctrl = 50)
  
  message(sprintf("    RawMean(NHAS) 与 nFeature_Spatial 相关: rho = %.3f (说明受 UMI 混杂极重!)",
                  cor(seu3$NHAS_RawMean, seu3$nFeature_Spatial, method = "spearman")))
  message(sprintf("    AddModuleScore(NHAS) 与 nFeature_Spatial 相关: rho = %.3f",
                  cor(seu3$NHAS_MS_1, seu3$nFeature_Spatial, method = "spearman")))
  message(sprintf("    去混杂后 AddModuleScore 相关 -> NHAS vs CTRL_G2M: %.3f | vs HM_Complement: %.3f | vs T_exhaustion: %.3f",
                  cor(seu3$NHAS_MS_1, seu3$SIG_MS_1, method = "spearman"),
                  cor(seu3$NHAS_MS_1, seu3$SIG_MS_2, method = "spearman"),
                  cor(seu3$NHAS_MS_1, seu3$SIG_MS_3, method = "spearman")))
  
  id_csv <- list.files(h3_dir, pattern = "identity\\.csv$", full.names = TRUE)[1]
  if (!is.na(id_csv)) {
    sme_df <- fread(id_csv, data.table = FALSE)
    message("    HCC3R 作者附带注释文件: ", basename(id_csv), " | 列名: ", paste(colnames(sme_df), collapse = ", "))
    print(head(sme_df, 3))
  }
}

message("\n>>> 全部完成！请查看 04_figures/Fig5_RealWorld_Clinical_Bias.png 并把控制台输出发给我。")
