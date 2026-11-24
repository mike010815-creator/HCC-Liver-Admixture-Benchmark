# ==============================================================================
# 脚本名称: 07b_interaction_clinical_spatial.R (修复完整版)
# 核心任务:
#   1. 深度排查基质交互项 (Pan_F_TBRs / HM_EMT / CAF_Peng) 的符号矛盾:
#      NHAS 偏度诊断 + 线性/Winsorize/秩正态化/二分交互对比 + Q1-Q4 四分位 HR 轨迹
#   2. 预设精简临床多变量 Cox 敏感性分析 (Age + Sex + TNM Stage III-IV vs I-II):
#      打印完整临床字典与缺失审计 + 固定/随机 Meta + C-index + PH 检验
#   3. GSE238264 全部 7 张 Visium 空间切片统一流程分析 (消除选片偏倚)
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
})

out_dir <- "03_results/revision_20261123"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tcga_obj <- readRDS("01_clean_data/tcga_lihc_clean.rds")
h18_obj  <- readRDS("01_clean_data/hccdb18_clean.rds")
gse_obj  <- readRDS("01_clean_data/gse14520_clean.rds")
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
bulk_sim <- readRDS("01_clean_data/02d_multicohort_sim_results.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
res4_cli <- readRDS("01_clean_data/04b_clinical_survival_results.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

sig_list     <- sigs_obj$sig_list
pur_sig_list <- res5_obj$purified_sig_list
pruned_sigs  <- res5_obj$prune_summary$Signature[res5_obj$prune_summary$Was_Pruned]
nhas_genes   <- res3_obj$nhas_genes

# 固定效应 + DL 随机效应 Meta 汇总函数
pool_meta <- function(beta, se) {
  ok <- !is.na(beta) & !is.na(se) & se > 0 & is.finite(beta) & is.finite(se)
  b <- beta[ok]; s <- se[ok]; k <- length(b)
  if (k == 0) {
    return(list(k = 0, HR = NA_real_, Lo = NA_real_, Hi = NA_real_, P = NA_real_,
                I2 = NA_real_, R_HR = NA_real_, R_Lo = NA_real_, R_Hi = NA_real_, R_P = NA_real_))
  }
  w <- 1 / (s^2); B <- sum(w * b) / sum(w); S <- sqrt(1 / sum(w))
  Q <- sum(w * (b - B)^2)
  I2 <- if (k > 1 && Q > 0) max(0, (Q - k + 1) / Q) * 100 else 0
  tau2 <- if (k > 1) max(0, (Q - k + 1) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (s^2 + tau2); Br <- sum(wr * b) / sum(wr); Sr <- sqrt(1 / sum(wr))
  list(
    k = k, HR = exp(B), Lo = exp(B - 1.96 * S), Hi = exp(B + 1.96 * S),
    P = 2 * pnorm(-abs(B / S)), I2 = I2,
    R_HR = exp(Br), R_Lo = exp(Br - 1.96 * Sr), R_Hi = exp(Br + 1.96 * Sr),
    R_P = 2 * pnorm(-abs(Br / Sr))
  )
}

# ==============================================================================
# 模块 4: 基质交互项“符号矛盾”深度排雷
# ==============================================================================
message(">>> [1/3] 正在排查基质交互项 (Pan_F_TBRs / HM_EMT / CAF_Peng) 的符号矛盾...")

nhas_dist <- map_dfr(names(res4_cli$real_all_list), function(cn) {
  v <- as.numeric(scale(res4_cli$real_all_list[[cn]]$NHAS))
  skew <- mean(v^3)
  qs <- quantile(v, c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
  tibble(Cohort = cn, Min = min(v), P01 = qs[1], P05 = qs[2], Median = qs[4],
         P95 = qs[6], Max = max(v), Skewness = skew)
})
message("    各队列 NHAS_z 分布与偏度 (负偏度代表存在极低分化左尾杠杆点):")
print(as.data.frame(nhas_dist %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)

diag_sigs <- c("Pan_F_TBRs", "HM_EMT", "CAF_Peng", "HM_Angiogenesis", "HM_TNFA_via_NFKB", "T_cell_inflamed_GEP")

int_diag_cohort <- map_dfr(names(res4_cli$real_all_list), function(cn) {
  dat_c <- res4_cli$real_all_list[[cn]]
  z_raw <- as.numeric(scale(dat_c$NHAS))
  q05 <- quantile(z_raw, 0.05); q95 <- quantile(z_raw, 0.95)
  dat_c$NHAS_raw <- z_raw
  dat_c$NHAS_win <- pmin(pmax(z_raw, q05), q95)
  dat_c$NHAS_rn  <- qnorm((rank(dat_c$NHAS) - 0.5) / nrow(dat_c))
  dat_c$NHAS_hi  <- as.numeric(dat_c$NHAS >= median(dat_c$NHAS))
  
  map_dfr(diag_sigs, function(s) {
    dat_s <- dat_c
    dat_s$Sig_z <- as.numeric(scale(dat_s[[s]]))
    sm1 <- summary(coxph(Surv(OS_time, OS_event) ~ Sig_z * NHAS_raw, data = dat_s))$coefficients
    sm2 <- summary(coxph(Surv(OS_time, OS_event) ~ Sig_z * NHAS_win, data = dat_s))$coefficients
    sm3 <- summary(coxph(Surv(OS_time, OS_event) ~ Sig_z * NHAS_rn,  data = dat_s))$coefficients
    sm4 <- summary(coxph(Surv(OS_time, OS_event) ~ Sig_z * NHAS_hi,  data = dat_s))$coefficients
    tibble(
      Cohort    = cn,
      Signature = s,
      Int_Type  = c("1_Linear_Raw_NHASz", "2_Winsor_5_95_NHASz", "3_RankNormal_NHAS", "4_Binary_High_vs_Low"),
      beta      = c(unname(sm1["Sig_z:NHAS_raw", "coef"]),
                    unname(sm2["Sig_z:NHAS_win", "coef"]),
                    unname(sm3["Sig_z:NHAS_rn",  "coef"]),
                    unname(sm4["Sig_z:NHAS_hi",  "coef"])),
      se        = c(unname(sm1["Sig_z:NHAS_raw", "se(coef)"]),
                    unname(sm2["Sig_z:NHAS_win", "se(coef)"]),
                    unname(sm3["Sig_z:NHAS_rn",  "se(coef)"]),
                    unname(sm4["Sig_z:NHAS_hi",  "se(coef)"]))
    )
  })
})

int_diag_meta <- int_diag_cohort %>%
  group_by(Signature, Int_Type) %>%
  summarise(res = list(pool_meta(beta, se)), .groups = "drop") %>%
  unnest_wider(res) %>%
  arrange(Signature, Int_Type)

message("\n    [交互模型敏感性诊断] 4 种交互定义对比:")
print(as.data.frame(
  int_diag_meta %>%
    mutate(across(c(HR, Lo, Hi, I2), ~ round(.x, 3)), P = signif(P, 3)) %>%
    select(Signature, Int_Type, HR, Lo, Hi, P, I2)
), row.names = FALSE)

# NHAS 四分位 (Q1 -> Q4) 分层 Cox (全队列统一 SD)
q_cohort <- map_dfr(names(res4_cli$real_all_list), function(cn) {
  dat_c <- res4_cli$real_all_list[[cn]]
  qs <- quantile(dat_c$NHAS, probs = c(0, 0.25, 0.50, 0.75, 1.0))
  dat_c$NHAS_Q <- cut(dat_c$NHAS, breaks = qs, include.lowest = TRUE,
                      labels = c("Q1_Lowest", "Q2", "Q3", "Q4_Highest"))
  map_dfr(c("Pan_F_TBRs", "HM_EMT", "CAF_Peng"), function(s) {
    dat_s <- dat_c
    dat_s$Sig_z <- as.numeric(scale(dat_s[[s]]))
    map_dfr(levels(dat_s$NHAS_Q), function(ql) {
      sub_d <- dat_s[dat_s$NHAS_Q == ql, , drop = FALSE]
      fit <- tryCatch(summary(coxph(Surv(OS_time, OS_event) ~ Sig_z, data = sub_d))$coefficients,
                      error = function(e) NULL)
      if (is.null(fit)) return(tibble(Cohort = cn, Signature = s, NHAS_Q = ql, beta = NA_real_, se = NA_real_))
      tibble(Cohort = cn, Signature = s, NHAS_Q = ql,
             beta = unname(fit[1, "coef"]), se = unname(fit[1, "se(coef)"]))
    })
  })
})

q_meta <- q_cohort %>%
  group_by(Signature, NHAS_Q) %>%
  summarise(res = list(pool_meta(beta, se)), .groups = "drop") %>%
  unnest_wider(res)

message("\n    [NHAS 四分位 HR 轨迹 (Q1 最低 -> Q4 最高, 统一全队列 SD)]:")
print(as.data.frame(
  q_meta %>% mutate(across(c(HR, Lo, Hi), ~ round(.x, 3)), P = signif(P, 3)) %>%
    select(Signature, NHAS_Q, HR, Lo, Hi, P)
), row.names = FALSE)

write.csv(int_diag_meta, file.path(out_dir, "07b_interaction_sign_diagnosis.csv"), row.names = FALSE)
write.csv(q_meta,        file.path(out_dir, "07b_stromal_quartile_HRs.csv"),       row.names = FALSE)

# ==============================================================================
# 模块 5: 预设精简临床多变量 Cox 敏感性分析 (Age + Sex + TNM Stage III-IV vs I-II)
# ==============================================================================
message("\n>>> [2/3] 正在提取三队列临床协变量并运行预设精简多变量 Cox (Age + Sex + Stage)...")

# 5.1 TCGA-LIHC
tcga_m <- tcga_obj$all_T_meta
col_t_age   <- grep("age_at_index|age_at_initial_pathologic_diagnosis|age_at_diagnosis",
                    colnames(tcga_m), ignore.case = TRUE, value = TRUE)[1]
col_t_sex   <- grep("^gender|gender\\.demographic|^sex$", colnames(tcga_m), ignore.case = TRUE, value = TRUE)[1]
col_t_stage <- grep("ajcc_pathologic_stage|tumor_stage", colnames(tcga_m), ignore.case = TRUE, value = TRUE)[1]
message(sprintf("    [TCGA 临床映射列] Age=%s | Sex=%s | Stage=%s", col_t_age, col_t_sex, col_t_stage))

t_age_num <- suppressWarnings(as.numeric(tcga_m[[col_t_age]]))
if (median(t_age_num, na.rm = TRUE) > 200) t_age_num <- t_age_num / 365.25

tcga_clin <- tibble(
  sample_id = tcga_m$sample_id,
  Age_raw   = t_age_num,
  Sex_raw   = trimws(as.character(tcga_m[[col_t_sex]])),
  Stage_raw = trimws(as.character(tcga_m[[col_t_stage]]))
) %>%
  mutate(
    Sex_Male   = case_when(grepl("^male$|^m$", Sex_raw, ignore.case = TRUE) ~ 1,
                           grepl("^female$|^f$", Sex_raw, ignore.case = TRUE) ~ 0, TRUE ~ NA_real_),
    Stage_High = case_when(grepl("III|IV", Stage_raw, ignore.case = TRUE) ~ 1,
                           grepl("Stage\\s*I($|[AB\\s])|Stage\\s*II($|[AB\\s])|^I($|[AB])|^II($|[AB])",
                                 Stage_raw, ignore.case = TRUE) ~ 0,
                           TRUE ~ NA_real_)
  ) %>% distinct(sample_id, .keep_all = TRUE)

# 5.2 ICGC-LIRI-JP (HCCDB18)
h18_t_samp <- h18_obj$samp_meta %>% filter(sample_id %in% colnames(h18_obj$all_T_tpm))
h18_pat    <- h18_obj$pat_meta
ov_id <- sum(h18_t_samp$pat_key %in% h18_pat$patient_id)
ov_do <- if ("PATIENT" %in% colnames(h18_t_samp) && "PATIENT" %in% colnames(h18_pat)) {
  sum(h18_t_samp$PATIENT %in% h18_pat$PATIENT)
} else 0
h18_m <- if (ov_id >= ov_do) {
  left_join(h18_t_samp, h18_pat, by = c("pat_key" = "patient_id"))
} else {
  left_join(h18_t_samp, h18_pat, by = "PATIENT")
}

message("    [ICGC 原始 Stage 取值]: TNM_STAGE_T -> ",
        paste(names(table(h18_m$TNM_STAGE_T)), table(h18_m$TNM_STAGE_T), sep = "=", collapse = " | "))
if ("TNM_STAGE_T1" %in% colnames(h18_m)) {
  message("                            TNM_STAGE_T1 -> ",
          paste(names(table(h18_m$TNM_STAGE_T1)), table(h18_m$TNM_STAGE_T1), sep = "=", collapse = " | "))
}

h18_clin <- tibble(
  sample_id = h18_m$sample_id,
  Age_raw   = suppressWarnings(as.numeric(h18_m$AGE)),
  Sex_raw   = trimws(as.character(h18_m$GENDER)),
  Stage_raw = trimws(as.character(h18_m$TNM_STAGE_T))
) %>%
  mutate(
    Sex_Male   = case_when(grepl("^male$|^m$", Sex_raw, ignore.case = TRUE) ~ 1,
                           grepl("^female$|^f$", Sex_raw, ignore.case = TRUE) ~ 0, TRUE ~ NA_real_),
    Stage_High = case_when(grepl("III|IV|^[34]$|T[34]", Stage_raw, ignore.case = TRUE) ~ 1,
                           grepl("^[12]$|T[12]|\\bI\\b|\\bII\\b|Stage\\s*I|Stage\\s*II", Stage_raw, ignore.case = TRUE) ~ 0,
                           TRUE ~ NA_real_)
  ) %>% distinct(sample_id, .keep_all = TRUE)

# 5.3 GSE14520
gse_m <- gse_obj$meta %>% filter(sample_id %in% colnames(gse_obj$all_T_lin))
col_g_stage <- grep("TNM", colnames(gse_m), ignore.case = TRUE, value = TRUE)[1]
message("    [GSE14520 原始 Stage 取值]: ",
        paste(names(table(gse_m[[col_g_stage]])), table(gse_m[[col_g_stage]]), sep = "=", collapse = " | "))

gse_clin <- tibble(
  sample_id = gse_m$sample_id,
  Age_raw   = suppressWarnings(as.numeric(gse_m$Age)),
  Sex_raw   = trimws(as.character(gse_m$Gender)),
  Stage_raw = trimws(as.character(gse_m[[col_g_stage]]))
) %>%
  mutate(
    Sex_Male   = case_when(grepl("^male$|^m$", Sex_raw, ignore.case = TRUE) ~ 1,
                           grepl("^female$|^f$", Sex_raw, ignore.case = TRUE) ~ 0, TRUE ~ NA_real_),
    Stage_High = case_when(grepl("III|IV|^[34]", Stage_raw, ignore.case = TRUE) ~ 1,
                           grepl("^I|^II|^[12]", Stage_raw, ignore.case = TRUE) ~ 0,
                           TRUE ~ NA_real_)
  ) %>% distinct(sample_id, .keep_all = TRUE)

clin_map <- list(TCGA_LIHC = tcga_clin, ICGC_LIRI_JP = h18_clin, GSE14520 = gse_clin)

clin_audit <- map_dfr(names(clin_map), function(cn) {
  d <- res4_cli$real_all_list[[cn]] %>% left_join(clin_map[[cn]], by = "sample_id")
  cc <- d %>% filter(!is.na(Age_raw), !is.na(Sex_Male), !is.na(Stage_High))
  tibble(
    Cohort = cn, N_Total = nrow(d), Events_Total = sum(d$OS_event),
    Miss_Age = sum(is.na(d$Age_raw)), Miss_Sex = sum(is.na(d$Sex_Male)), Miss_Stage = sum(is.na(d$Stage_High)),
    N_Complete = nrow(cc), Events_Complete = sum(cc$OS_event),
    Pct_Stage_High = round(mean(cc$Stage_High) * 100, 1),
    EPV_4df = round(sum(cc$OS_event) / 4, 1)
  )
})
message("\n    [三队列临床协变量完整度与 EPV 审计表]:")
print(as.data.frame(clin_audit), row.names = FALSE)
write.csv(clin_audit, file.path(out_dir, "07b_clinical_covariate_audit.csv"), row.names = FALSE)

# 5.4 多变量 Cox (Age + Sex + Stage) 与 C-index
run_ss_mat <- function(lin_mat, gsets) {
  gsva(ssgseaParam(as.matrix(log2(pmax(lin_mat, 0) + 1)), gsets, minSize = 2, normalize = FALSE), verbose = FALSE)
}
real_mats <- list(TCGA_LIHC = tcga_obj$all_T_tpm, ICGC_LIRI_JP = h18_obj$all_T_tpm, GSE14520 = gse_obj$all_T_lin)
key_pruned_check <- c("HM_Complement", "HM_IL2_STAT5", "T_cell_exhaustion",
                      "MDSC_Peng", "HM_IFN_alpha", "HM_IFN_gamma", "TAM_Peng")

mv_pruned_cohort <- map_dfr(names(clin_map), function(cn) {
  d <- res4_cli$real_all_list[[cn]] %>%
    left_join(clin_map[[cn]], by = "sample_id") %>%
    filter(!is.na(Age_raw), !is.na(Sex_Male), !is.na(Stage_High))
  d$Age_z  <- as.numeric(scale(d$Age_raw))
  d$NHAS_z <- as.numeric(scale(d$NHAS))
  
  m_cc <- real_mats[[cn]][, d$sample_id, drop = FALSE]
  sc_p <- run_ss_mat(m_cc, pur_sig_list[pruned_sigs])
  
  f_clin <- coxph(Surv(OS_time, OS_event) ~ Age_z + Sex_Male + Stage_High, data = d)
  c_clin <- unname(summary(f_clin)$concordance[1])
  
  map_dfr(pruned_sigs, function(s) {
    ds <- d
    ds$Raw_z <- as.numeric(scale(ds[[s]]))
    ds$Pur_z <- as.numeric(scale(as.numeric(sc_p[s, ])))
    
    f_uni_p  <- coxph(Surv(OS_time, OS_event) ~ Pur_z, data = ds)
    f_mv_p   <- coxph(Surv(OS_time, OS_event) ~ Pur_z + Age_z + Sex_Male + Stage_High, data = ds)
    f_mv_r   <- coxph(Surv(OS_time, OS_event) ~ Raw_z + Age_z + Sex_Male + Stage_High, data = ds)
    f_mv_pnh <- coxph(Surv(OS_time, OS_event) ~ Pur_z + NHAS_z + Age_z + Sex_Male + Stage_High, data = ds)
    
    sp_u  <- summary(f_uni_p)$coefficients["Pur_z", ]
    sp_mv <- summary(f_mv_p)$coefficients["Pur_z", ]
    sr_mv <- summary(f_mv_r)$coefficients["Raw_z", ]
    sp_nh <- summary(f_mv_pnh)$coefficients["Pur_z", ]
    ph_p  <- tryCatch(cox.zph(f_mv_p)$table["GLOBAL", "p"], error = function(e) NA_real_)
    
    tibble(
      Cohort = cn, Signature = s,
      beta_pur_uni = unname(sp_u["coef"]),   se_pur_uni = unname(sp_u["se(coef)"]),
      beta_pur_mv  = unname(sp_mv["coef"]),  se_pur_mv  = unname(sp_mv["se(coef)"]),
      beta_raw_mv  = unname(sr_mv["coef"]),  se_raw_mv  = unname(sr_mv["se(coef)"]),
      beta_pur_mv_nhas = unname(sp_nh["coef"]), se_pur_mv_nhas = unname(sp_nh["se(coef)"]),
      Cidx_Clin     = c_clin,
      Cidx_Clin_Raw = unname(summary(f_mv_r)$concordance[1]),
      Cidx_Clin_Pur = unname(summary(f_mv_p)$concordance[1]),
      Global_PH_P   = unname(ph_p)
    )
  })
})

mv_pruned_meta <- mv_pruned_cohort %>%
  group_by(Signature) %>%
  summarise(
    u_res  = list(pool_meta(beta_pur_uni, se_pur_uni)),
    mv_res = list(pool_meta(beta_pur_mv,  se_pur_mv)),
    rm_res = list(pool_meta(beta_raw_mv,  se_raw_mv)),
    Mean_dC_Raw = mean(Cidx_Clin_Raw - Cidx_Clin),
    Mean_dC_Pur = mean(Cidx_Clin_Pur - Cidx_Clin),
    .groups = "drop"
  ) %>%
  mutate(
    Raw_MV_HR      = sapply(rm_res, `[[`, "HR"),   Raw_MV_P      = sapply(rm_res, `[[`, "P"),
    Pur_Uni_HR     = sapply(u_res,  `[[`, "HR"),   Pur_Uni_P     = sapply(u_res,  `[[`, "P"),
    Pur_MV_HR      = sapply(mv_res, `[[`, "HR"),   Pur_MV_Lo     = sapply(mv_res, `[[`, "Lo"),
    Pur_MV_Hi      = sapply(mv_res, `[[`, "Hi"),   Pur_MV_P      = sapply(mv_res, `[[`, "P"),
    Pur_MV_I2      = sapply(mv_res, `[[`, "I2"),   Pur_MV_Rand_HR = sapply(mv_res, `[[`, "R_HR"),
    Pur_MV_Rand_P  = sapply(mv_res, `[[`, "R_P")
  ) %>%
  select(-ends_with("_res")) %>%
  mutate(Pur_MV_BH_q = p.adjust(Pur_MV_P, "BH")) %>%
  arrange(Pur_MV_P)

write.csv(mv_pruned_cohort, file.path(out_dir, "07b_multivar_pruned_cohort.csv"), row.names = FALSE)
write.csv(mv_pruned_meta,   file.path(out_dir, "07b_multivar_pruned_meta.csv"),   row.names = FALSE)

message("\n    [结果 2] 调整 Age + Sex + TNM Stage 后的多变量 Meta-Cox (S3 净化 vs 原始):")
print(as.data.frame(
  mv_pruned_meta %>%
    filter(Signature %in% key_pruned_check | Pur_MV_P < 0.05) %>%
    mutate(across(c(Raw_MV_HR, Pur_MV_HR, Pur_MV_Lo, Pur_MV_Hi, Pur_MV_Rand_HR), ~ round(.x, 3)),
           across(c(Mean_dC_Raw, Mean_dC_Pur), ~ round(.x, 4)),
           Pur_MV_I2 = round(Pur_MV_I2, 1),
           across(c(Raw_MV_P, Pur_MV_P, Pur_MV_BH_q, Pur_MV_Rand_P), ~ signif(.x, 3))) %>%
    select(Signature, Raw_MV_HR, Raw_MV_P, Pur_MV_HR, Pur_MV_Lo, Pur_MV_Hi,
           Pur_MV_P, Pur_MV_BH_q, Pur_MV_I2, Pur_MV_Rand_P, Mean_dC_Pur)
), row.names = FALSE)

# ==============================================================================
# 模块 6: GSE238264 全部 7 张 Visium 空间切片统一分析 (消除选片偏倚)
# ==============================================================================
message("\n>>> [3/3] 正在对 GSE238264 全部 7 张 Visium 空间切片执行统一 AddModuleScore 流程...")

sp_root <- "01_clean_data/spatial_unpacked/GSE238264"
slide_dirs <- list.dirs(sp_root, recursive = FALSE, full.names = TRUE)
slide_dirs <- slide_dirs[file.exists(file.path(slide_dirs, "spatial", "scalefactors_json.json"))]
slide_names <- sort(basename(slide_dirs))

score_ms_safe <- function(obj, gvec, tag) {
  gu <- intersect(gvec, rownames(obj))
  if (length(gu) < 2) return(rep(NA_real_, ncol(obj)))
  tmp <- AddModuleScore(obj, features = list(gu), name = tag, ctrl = 50, seed = 2026)
  as.numeric(tmp[[paste0(tag, "1")]][, 1])
}

bulk_mean30 <- bulk_sim %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  group_by(Signature) %>%
  summarise(Bulk30_Mean = mean(Std_Drift), .groups = "drop")

sp7_qc_list  <- list()
sp7_cor_list <- list()

for (sn in slide_names) {
  s_dir <- file.path(sp_root, sn)
  cnt <- Read10X(file.path(s_dir, "filtered_feature_bc_matrix"))
  if (is.list(cnt)) cnt <- cnt[[1]]
  seu <- CreateSeuratObject(counts = cnt, assay = "Spatial", project = sn)
  img <- Read10X_Image(file.path(s_dir, "spatial"), image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
  DefaultAssay(img) <- "Spatial"
  bc <- intersect(colnames(seu), Cells(img))
  seu <- subset(seu, cells = bc)
  seu <- subset(seu, subset = nFeature_Spatial >= 500)
  seu <- NormalizeData(seu, verbose = FALSE)
  
  nhas_ms <- score_ms_safe(seu, nhas_genes, "M_NHAS_")
  raw_ms  <- sapply(names(sig_list), function(s) score_ms_safe(seu, sig_list[[s]], "M_RAW_"))
  pur_ms  <- sapply(pruned_sigs,     function(s) score_ms_safe(seu, pur_sig_list[[s]], "M_PUR_"))
  
  raw_rhos <- apply(raw_ms, 2, function(v) suppressWarnings(cor(nhas_ms, v, method = "spearman", use = "complete.obs")))
  pur_rhos <- apply(pur_ms, 2, function(v) suppressWarnings(cor(nhas_ms, v, method = "spearman", use = "complete.obs")))
  
  b_aligned <- bulk_mean30$Bulk30_Mean[match(names(raw_rhos), bulk_mean30$Signature)]
  conc_bulk <- cor(raw_rhos, b_aligned, method = "spearman", use = "complete.obs")
  q_n <- quantile(nhas_ms, c(0.10, 0.90), na.rm = TRUE)
  
  sp7_qc_list[[sn]] <- tibble(
    Slide = sn,
    N_Spots_ge500 = ncol(seu),
    NHAS_MS_P10 = round(q_n[1], 3),
    NHAS_MS_P90 = round(q_n[2], 3),
    NHAS_MS_Span = round(q_n[2] - q_n[1], 3),
    Residual_Rho_nFeature = round(cor(nhas_ms, seu$nFeature_Spatial, method = "spearman"), 3),
    Rho_NHAS_vs_G2M = round(raw_rhos["CTRL_G2M_Checkpoint"], 3),
    Rho_NHAS_vs_Complement = round(raw_rhos["HM_Complement"], 3),
    Concordance_vs_Bulk30 = round(conc_bulk, 3)
  )
  
  sp7_cor_list[[sn]] <- tibble(
    Slide = sn,
    Signature = names(raw_rhos),
    Spatial_Rho_Raw = as.numeric(raw_rhos),
    Spatial_Rho_Pur = as.numeric(pur_rhos[match(names(raw_rhos), names(pur_rhos))])
  )
  message(sprintf("    %s: %d spots | NHAS span=%.2f | 残余深度rho=%.2f | vs Bulk rho=%.3f",
                  sn, ncol(seu), q_n[2] - q_n[1],
                  cor(nhas_ms, seu$nFeature_Spatial, method = "spearman"), conc_bulk))
}

sp7_qc_df  <- bind_rows(sp7_qc_list) %>% arrange(desc(NHAS_MS_Span))
sp7_cor_df <- bind_rows(sp7_cor_list)

sp7_sig_summary <- sp7_cor_df %>%
  group_by(Signature) %>%
  summarise(
    Median_Rho_7slides = median(Spatial_Rho_Raw, na.rm = TRUE),
    IQR_Rho_7slides    = IQR(Spatial_Rho_Raw, na.rm = TRUE),
    N_Slides_Positive  = sum(Spatial_Rho_Raw > 0, na.rm = TRUE),
    Median_Pur_7slides = median(Spatial_Rho_Pur, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(bulk_mean30, by = "Signature") %>%
  arrange(desc(Median_Rho_7slides))

conc_7slide_med <- cor(sp7_sig_summary$Median_Rho_7slides, sp7_sig_summary$Bulk30_Mean, method = "spearman")

write.csv(sp7_qc_df,       file.path(out_dir, "07b_spatial_all7_slides_qc.csv"),  row.names = FALSE)
write.csv(sp7_cor_df,      file.path(out_dir, "07b_spatial_all7_slides_cor.csv"), row.names = FALSE)
write.csv(sp7_sig_summary, file.path(out_dir, "07b_spatial_all7_sig_summary.csv"), row.names = FALSE)

message("\n>>> [结果 3A] 全部 7 张 Visium 空间切片统一质控与 Bulk 一致性总表:")
print(as.data.frame(sp7_qc_df), row.names = FALSE)

message(sprintf("\n>>> [结果 3B] 7 切片中位空间 NHAS 相关 vs Bulk 30%% 漂移总体一致性: Spearman rho = %.3f",
                conc_7slide_med))
print(as.data.frame(
  sp7_sig_summary %>% mutate(across(c(Median_Rho_7slides, IQR_Rho_7slides, Median_Pur_7slides, Bulk30_Mean), ~ round(.x, 3)))
), row.names = FALSE)

message("\n>>> 07b 全部完成！所有补算表格已保存至 ", out_dir)
