# ==============================================================================
# 脚本名称: 07a_bootstrap_null_and_external_S2.R
# 核心任务:
#   1. 患者配对 Bootstrap (B = 1000): 给出 Raw vs S3 在自身基线与共同原始基线上的
#      精确 95% 置信区间，并按基线保真度 (>=0.90 vs <0.90) 分层剖析
#   2. 严格零 Bulk 泄漏的随机删基因阴性对照 (R = 50 均匀随机 + R = 50 表达分箱匹配)
#   3. 严格跨队列外部参考 S2 验证 (ICGC -> TCGA; TCGA -> ICGC)
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
})

out_dir <- "03_results/revision_20261123"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. 加载已清洗数据与第 3/5 关结果
# ------------------------------------------------------------------------------
tcga_obj <- readRDS("01_clean_data/tcga_lihc_clean.rds")
h18_obj  <- readRDS("01_clean_data/hccdb18_clean.rds")
gse_obj  <- readRDS("01_clean_data/gse14520_clean.rds")
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

sig_list      <- sigs_obj$sig_list
pur_sig_list  <- res5_obj$purified_sig_list
prune_summary <- res5_obj$prune_summary
pruned_sigs   <- prune_summary$Signature[prune_summary$Was_Pruned]
nhas_genes    <- res3_obj$nhas_genes
sig_gene_map  <- res3_obj$sig_gene_map

cohorts <- list(
  TCGA_LIHC    = list(T_mat = tcga_obj$T_tpm, N_mat = tcga_obj$N_tpm),
  ICGC_LIRI_JP = list(T_mat = h18_obj$T_tpm,  N_mat = h18_obj$N_tpm),
  GSE14520     = list(T_mat = gse_obj$T_lin,  N_mat = gse_obj$N_lin)
)

run_ss <- function(lin_mat, gsets) {
  m <- log2(pmax(lin_mat, 0) + 1)
  gsva(ssgseaParam(as.matrix(m), gsets, minSize = 2, normalize = FALSE), verbose = FALSE)
}

# 高速精确逐行 Spearman 相关 (避免循环调用 cor.test，误差 < 1e-15)
row_spearman_fast <- function(A, B) {
  n <- ncol(A)
  rA <- apply(t(A), 2, rank, ties.method = "average")
  rB <- apply(t(B), 2, rank, ties.method = "average")
  sA <- scale(rA)
  sB <- scale(rB)
  vals <- colSums(sA * sB) / (n - 1)
  as.numeric(vals)
}

# ==============================================================================
# 模块 1: 患者配对 Bootstrap (B = 1000) —— 自身基线 & 共同原始目标 95% CI
# ==============================================================================
message(">>> [1/3] 正在计算三队列 30% 混入单样本得分并运行患者配对 Bootstrap (B = 1000)...")

scores_30 <- list()
for (cn in names(cohorts)) {
  T_mat <- cohorts[[cn]]$T_mat
  N_mat <- cohorts[[cn]]$N_mat
  mix30 <- 0.70 * T_mat + 0.30 * N_mat
  scores_30[[cn]] <- list(
    r0  = run_ss(T_mat, sig_list[pruned_sigs]),
    r30 = run_ss(mix30, sig_list[pruned_sigs]),
    p0  = run_ss(T_mat, pur_sig_list[pruned_sigs]),
    p30 = run_ss(mix30, pur_sig_list[pruned_sigs])
  )
}

# 划分高基线一致 (Fidelity >= 0.90) vs 重度重构 (< 0.90) 两个特征子组
fid_df <- res5_obj$fidelity_df %>% filter(Signature %in% pruned_sigs)
high_fid_sigs <- fid_df$Signature[fid_df$Mean_PureTumor_Fidelity_Rho >= 0.90]
low_fid_sigs  <- fid_df$Signature[fid_df$Mean_PureTumor_Fidelity_Rho <  0.90]
message(sprintf("    高基线一致特征 (Fidelity >= 0.90, n=%d): %s",
                length(high_fid_sigs), paste(high_fid_sigs, collapse = ", ")))
message(sprintf("    重度重构特征   (Fidelity <  0.90, n=%d): %s",
                length(low_fid_sigs),  paste(low_fid_sigs,  collapse = ", ")))

calc_cohort_vec <- function(sc, idx) {
  r0  <- sc$r0[,  idx, drop = FALSE]
  r30 <- sc$r30[, idx, drop = FALSE]
  p0  <- sc$p0[,  idx, drop = FALSE]
  p30 <- sc$p30[, idx, drop = FALSE]
  
  sd_r0 <- apply(r0, 1, sd); sd_p0 <- apply(p0, 1, sd)
  med_r0  <- apply(r0,  1, median); med_r30 <- apply(r30, 1, median)
  med_p0  <- apply(p0,  1, median); med_p30 <- apply(p30, 1, median)
  
  list(
    raw_abs_drift   = abs(rowMeans(r30 - r0) / sd_r0),
    s3_abs_drift    = abs(rowMeans(p30 - p0) / sd_p0),
    raw_rho         = row_spearman_fast(r0, r30),
    s3_self_rho     = row_spearman_fast(p0, p30),
    s3_orig_rho     = row_spearman_fast(r0, p30),
    raw_fflip       = rowMeans((r0 >= med_r0) != (r30 >= med_r0)) * 100,
    s3_self_fflip   = rowMeans((p0 >= med_p0) != (p30 >= med_p0)) * 100,
    raw_rflip       = rowMeans((r0 >= med_r0) != (r30 >= med_r30)) * 100,
    s3_self_rflip   = rowMeans((p0 >= med_p0) != (p30 >= med_p30)) * 100,
    s3_orig_rflip   = rowMeans((r0 >= med_r0) != (p30 >= med_p30)) * 100
  )
}

# 点估计 + B = 1000 患者配对重抽样
B <- 1000L
sig_subsets <- list(
  All_18_Pruned     = seq_along(pruned_sigs),
  High_Fidelity_ge09 = which(pruned_sigs %in% high_fid_sigs),
  Redesigned_lt09    = which(pruned_sigs %in% low_fid_sigs)
)

metric_names <- c(
  "raw_abs_drift", "s3_abs_drift", "delta_abs_drift",
  "raw_rho", "s3_self_rho", "delta_self_rho",
  "s3_orig_rho", "delta_orig_rho",
  "raw_fflip", "s3_self_fflip", "delta_fflip",
  "raw_rflip", "s3_self_rflip", "delta_self_rflip",
  "s3_orig_rflip", "delta_orig_rflip"
)

extract_summary_metrics <- function(cohort_res_list, s_idx) {
  # 等权三队列均值
  avg_sig <- function(k) {
    mat <- sapply(cohort_res_list, function(x) x[[k]][s_idx])
    mean(mat, na.rm = TRUE)
  }
  ra <- avg_sig("raw_abs_drift"); sa <- avg_sig("s3_abs_drift")
  rr <- avg_sig("raw_rho");       sr <- avg_sig("s3_self_rho"); sor <- avg_sig("s3_orig_rho")
  rf <- avg_sig("raw_fflip");     sf <- avg_sig("s3_self_fflip")
  rl <- avg_sig("raw_rflip");     sl <- avg_sig("s3_self_rflip"); sol <- avg_sig("s3_orig_rflip")
  c(
    raw_abs_drift = ra, s3_abs_drift = sa, delta_abs_drift = sa - ra,
    raw_rho = rr, s3_self_rho = sr, delta_self_rho = sr - rr,
    s3_orig_rho = sor, delta_orig_rho = sor - rr,
    raw_fflip = rf, s3_self_fflip = sf, delta_fflip = sf - rf,
    raw_rflip = rl, s3_self_rflip = sl, delta_self_rflip = sl - rl,
    s3_orig_rflip = sol, delta_orig_rflip = sol - rl
  )
}

# 点估计
obs_by_cohort <- lapply(names(cohorts), function(cn) {
  calc_cohort_vec(scores_30[[cn]], seq_len(ncol(scores_30[[cn]]$r0)))
})
names(obs_by_cohort) <- names(cohorts)

# Bootstrap 存储阵列
boot_store <- lapply(sig_subsets, function(x) matrix(NA_real_, nrow = B, ncol = length(metric_names),
                                                     dimnames = list(NULL, metric_names)))
# 同时保存按队列分开的 All_18 bootstrap
boot_by_cohort <- lapply(names(cohorts), function(cn) {
  matrix(NA_real_, nrow = B, ncol = length(metric_names), dimnames = list(NULL, metric_names))
})
names(boot_by_cohort) <- names(cohorts)

set.seed(2026)
for (b in seq_len(B)) {
  b_list <- lapply(names(cohorts), function(cn) {
    n_c <- ncol(scores_30[[cn]]$r0)
    idx <- sample.int(n_c, n_c, replace = TRUE)
    calc_cohort_vec(scores_30[[cn]], idx)
  })
  names(b_list) <- names(cohorts)
  
  for (sub_n in names(sig_subsets)) {
    boot_store[[sub_n]][b, ] <- extract_summary_metrics(b_list, sig_subsets[[sub_n]])
  }
  for (cn in names(cohorts)) {
    boot_by_cohort[[cn]][b, ] <- extract_summary_metrics(b_list[cn], sig_subsets$All_18_Pruned)
  }
}

summarize_boot <- function(obs_vec, bmat, group_lbl) {
  tibble(
    Group    = group_lbl,
    Metric   = names(obs_vec),
    Estimate = as.numeric(obs_vec),
    CI_Low   = apply(bmat, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_High  = apply(bmat, 2, quantile, probs = 0.975, na.rm = TRUE),
    Boot_P_2sided = sapply(seq_along(obs_vec), function(j) {
      v <- bmat[, j]
      2 * min(mean(v <= 0, na.rm = TRUE), mean(v >= 0, na.rm = TRUE))
    })
  )
}

boot_summary_df <- bind_rows(
  lapply(names(sig_subsets), function(sn) {
    summarize_boot(extract_summary_metrics(obs_by_cohort, sig_subsets[[sn]]), boot_store[[sn]], sn)
  }),
  lapply(names(cohorts), function(cn) {
    summarize_boot(extract_summary_metrics(obs_by_cohort[cn], sig_subsets$All_18_Pruned),
                   boot_by_cohort[[cn]], paste0("Cohort_", cn))
  })
)

write.csv(boot_summary_df, file.path(out_dir, "07a_patient_bootstrap_95CI.csv"), row.names = FALSE)

message("\n>>> [结果 1] 患者配对 Bootstrap 95% CI (均匀 30% 混入, B=1000, 全体 vs 高保真 vs 重构特征):")
print(as.data.frame(
  boot_summary_df %>%
    filter(Group %in% c("All_18_Pruned", "High_Fidelity_ge09", "Redesigned_lt09"),
           Metric %in% c("raw_abs_drift", "s3_abs_drift", "delta_abs_drift",
                         "raw_rho", "s3_self_rho", "delta_self_rho",
                         "s3_orig_rho", "delta_orig_rho",
                         "raw_fflip", "s3_self_fflip", "delta_fflip",
                         "raw_rflip", "s3_orig_rflip", "delta_orig_rflip")) %>%
    mutate(across(c(Estimate, CI_Low, CI_High), ~ round(.x, 3)),
           Boot_P_2sided = signif(Boot_P_2sided, 3))
), row.names = FALSE)

# ==============================================================================
# 模块 2: 匹配随机删基因阴性对照 (R = 50 均匀随机 + R = 50 表达量分箱匹配)
# ==============================================================================
message("\n>>> [2/3] 正在构造规模+表达量双匹配随机删基因对照 (R = 50/类, 零 Bulk 泄漏)...")

safe_sample <- function(x, size) {
  if (size <= 0L) return(character(0))
  if (length(x) == 1L) return(rep(x, size))
  sample(x, size, replace = FALSE)
}

R_null <- 50L
null_gsets_unif <- list()
null_gsets_expr <- list()

set.seed(2026)
for (s in pruned_sigs) {
  orig_g <- unique(sig_list[[s]])
  kept_g <- unique(pur_sig_list[[s]])
  K <- length(kept_g)
  
  # 用独立单细胞 GSE149614 原发肿瘤 Pseudo-bulk 表达量做分箱 (完全不接触 Bulk 数据)
  sub_sc <- sig_gene_map %>% filter(Signature == s)
  cpm_vec <- setNames(rep(0, length(orig_g)), orig_g)
  hit_g <- intersect(orig_g, sub_sc$Gene)
  cpm_vec[hit_g] <- sub_sc$PB_Tumor_CPM[match(hit_g, sub_sc$Gene)]
  
  n_bins <- min(4L, max(2L, floor(length(orig_g) / 5L)))
  rk <- rank(cpm_vec, ties.method = "first")
  bins <- cut(rk, breaks = n_bins, labels = FALSE)
  bin_genes <- split(orig_g, bins)
  k_per_bin <- sapply(bin_genes, function(bg) sum(kept_g %in% bg))
  stopifnot(sum(k_per_bin) == K)
  
  for (r in seq_len(R_null)) {
    null_gsets_unif[[paste0(s, "__", r)]] <- safe_sample(orig_g, K)
    expr_pick <- unlist(mapply(safe_sample, bin_genes, k_per_bin, SIMPLIFY = FALSE), use.names = FALSE)
    null_gsets_expr[[paste0(s, "__", r)]] <- expr_pick
  }
}

eval_null_sets <- function(gset_list, null_type_lbl) {
  map_dfr(names(cohorts), function(cn) {
    T_mat <- cohorts[[cn]]$T_mat
    N_mat <- cohorts[[cn]]$N_mat
    mix30 <- 0.70 * T_mat + 0.30 * N_mat
    sc0  <- run_ss(T_mat, gset_list)
    sc30 <- run_ss(mix30, gset_list)
    r0_full <- scores_30[[cn]]$r0
    
    nm <- rownames(sc0)
    sig_vec <- sub("__[0-9]+$", "", nm)
    rep_vec <- as.integer(sub("^.*__", "", nm))
    
    sd0   <- apply(sc0, 1, sd)
    med0  <- apply(sc0, 1, median)
    med30 <- apply(sc30, 1, median)
    abs_d <- abs(rowMeans(sc30 - sc0) / sd0)
    self_rho <- row_spearman_fast(sc0, sc30)
    self_fflip <- rowMeans((sc0 >= med0) != (sc30 >= med0)) * 100
    
    # 与对应原始基因集基线的排序相关性
    r0_aligned <- r0_full[sig_vec, , drop = FALSE]
    orig_rho <- row_spearman_fast(r0_aligned, sc30)
    base_fid <- row_spearman_fast(r0_aligned, sc0)
    
    message(sprintf("    [%s] %s 计算完成 (%d 随机基因集)", null_type_lbl, cn, nrow(sc0)))
    tibble(
      Null_Type  = null_type_lbl,
      Cohort     = cn,
      Signature  = sig_vec,
      Rep        = rep_vec,
      Abs_Drift  = abs_d,
      Self_Rho   = self_rho,
      Self_FFlip = self_fflip,
      Orig_Rho   = orig_rho,
      Base_Fid   = base_fid
    )
  })
}

null_res_all <- bind_rows(
  eval_null_sets(null_gsets_unif, "Uniform_SizeMatched"),
  eval_null_sets(null_gsets_expr, "ExprBin_SizeMatched")
)

# 计算 S3 实际值 vs 50 次随机删基因零分布 (队列+18特征均值层面 & 逐特征层面)
s3_obs_by_sig <- map_dfr(names(cohorts), function(cn) {
  x <- obs_by_cohort[[cn]]
  tibble(
    Cohort    = cn,
    Signature = pruned_sigs,
    S3_Abs_Drift  = x$s3_abs_drift,
    S3_Self_Rho   = x$s3_self_rho,
    S3_Self_FFlip = x$s3_self_fflip,
    S3_Orig_Rho   = x$s3_orig_rho,
    S3_Base_Fid   = row_spearman_fast(scores_30[[cn]]$r0, scores_30[[cn]]$p0)
  )
})

s3_grand <- s3_obs_by_sig %>%
  summarise(across(starts_with("S3_"), mean))

null_rep_grand <- null_res_all %>%
  group_by(Null_Type, Rep) %>%
  summarise(
    Abs_Drift  = mean(Abs_Drift),
    Self_Rho   = mean(Self_Rho),
    Self_FFlip = mean(Self_FFlip),
    Orig_Rho   = mean(Orig_Rho),
    Base_Fid   = mean(Base_Fid),
    .groups    = "drop"
  )

null_grand_summary <- null_rep_grand %>%
  group_by(Null_Type) %>%
  summarise(
    S3_AbsDrift    = s3_grand$S3_Abs_Drift,
    Null_AbsDrift_Mean = mean(Abs_Drift),
    Null_AbsDrift_P05  = quantile(Abs_Drift, 0.05),
    Null_AbsDrift_P95  = quantile(Abs_Drift, 0.95),
    Emp_P_Drift_Lower  = (sum(Abs_Drift <= s3_grand$S3_Abs_Drift) + 1) / (n() + 1),
    
    S3_SelfFlip    = s3_grand$S3_Self_FFlip,
    Null_Flip_Mean = mean(Self_FFlip),
    Emp_P_Flip_Lower = (sum(Self_FFlip <= s3_grand$S3_Self_FFlip) + 1) / (n() + 1),
    
    S3_SelfRho     = s3_grand$S3_Self_Rho,
    Null_SelfRho_Mean = mean(Self_Rho),
    Emp_P_SelfRho_Higher = (sum(Self_Rho >= s3_grand$S3_Self_Rho) + 1) / (n() + 1),
    
    S3_OrigRho     = s3_grand$S3_Orig_Rho,
    Null_OrigRho_Mean = mean(Orig_Rho),
    S3_BaseFid     = s3_grand$S3_Base_Fid,
    Null_BaseFid_Mean = mean(Base_Fid),
    .groups = "drop"
  )

# 逐特征检验 (在表达量匹配对照中，有多少个特征 S3 显著优于随机删基因)
null_sig_summary <- null_res_all %>%
  filter(Null_Type == "ExprBin_SizeMatched") %>%
  group_by(Signature, Rep) %>%
  summarise(Abs_Drift = mean(Abs_Drift), Self_FFlip = mean(Self_FFlip),
            Self_Rho = mean(Self_Rho), Orig_Rho = mean(Orig_Rho), .groups = "drop") %>%
  left_join(
    s3_obs_by_sig %>% group_by(Signature) %>%
      summarise(across(starts_with("S3_"), mean), .groups = "drop"),
    by = "Signature"
  ) %>%
  group_by(Signature) %>%
  summarise(
    S3_Drift     = first(S3_Abs_Drift),
    Null_Drift   = mean(Abs_Drift),
    Emp_P_Drift  = (sum(Abs_Drift <= first(S3_Abs_Drift)) + 1) / (n() + 1),
    S3_Flip      = first(S3_Self_FFlip),
    Null_Flip    = mean(Self_FFlip),
    Emp_P_Flip   = (sum(Self_FFlip <= first(S3_Self_FFlip)) + 1) / (n() + 1),
    S3_OrigRho   = first(S3_Orig_Rho),
    Null_OrigRho = mean(Orig_Rho),
    .groups = "drop"
  ) %>%
  arrange(Emp_P_Drift, S3_Drift - Null_Drift)

write.csv(null_grand_summary, file.path(out_dir, "07a_random_deletion_null_grand.csv"), row.names = FALSE)
write.csv(null_sig_summary,   file.path(out_dir, "07a_random_deletion_null_bysig.csv"), row.names = FALSE)

message("\n>>> [结果 2A] S3 vs 随机删基因阴性对照总体对比 (18 特征 x 3 队列均值, R=50):")
print(as.data.frame(null_grand_summary %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)

message("\n>>> [结果 2B] 逐特征对比 (S3 vs 表达量+规模双匹配随机删基因):")
print(as.data.frame(null_sig_summary %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# ==============================================================================
# 模块 3: 严格跨队列外部参考 S2 验证 (ICGC -> TCGA; TCGA -> ICGC)
# ==============================================================================
message("\n>>> [3/3] 正在运行严格跨队列外部参考 S2 验证 (零训练泄漏)...")

purify_ext <- function(mix_lin, T_ref, N_ref, nhas_g) {
  common_g <- Reduce(intersect, list(rownames(mix_lin), rownames(T_ref), rownames(N_ref)))
  mix_c <- mix_lin[common_g, , drop = FALSE]
  mu_T  <- apply(T_ref[common_g, , drop = FALSE], 1, median)
  mu_N  <- apply(N_ref[common_g, , drop = FALSE], 1, median)
  
  g_use <- intersect(nhas_g, common_g)
  g_use <- g_use[mu_N[g_use] > 1.5 * mu_T[g_use] & mu_N[g_use] > 10]
  denom <- mu_N[g_use] - mu_T[g_use]
  
  est_mat <- sweep(sweep(mix_c[g_use, , drop = FALSE], 1, mu_T[g_use], "-"), 1, denom, "/")
  alpha_hat <- pmin(pmax(apply(est_mat, 2, median, na.rm = TRUE), 0), 0.65)
  
  pur_mat <- matrix(0, nrow = nrow(mix_c), ncol = ncol(mix_c), dimnames = dimnames(mix_c))
  for (j in seq_len(ncol(mix_c))) {
    aj <- alpha_hat[j]
    v <- (mix_c[, j] - aj * mu_N) / max(1 - aj, 0.35)
    v[v < 0] <- 0
    s_v <- sum(v)
    pur_mat[, j] <- if (s_v > 0) v / s_v * 1e6 else mix_c[, j]
  }
  list(pur_mat = pur_mat, alpha_hat = alpha_hat)
}

run_ext_s2_pair <- function(train_name, test_name) {
  T_tr <- cohorts[[train_name]]$T_mat; N_tr <- cohorts[[train_name]]$N_mat
  T_te <- cohorts[[test_name]]$T_mat;  N_te <- cohorts[[test_name]]$N_mat
  mix30 <- 0.70 * T_te + 0.30 * N_te
  
  # 同一公共基因全集上评估，防止基因宇宙差异带来偏差
  common_g <- Reduce(intersect, list(rownames(T_te), rownames(T_tr), rownames(N_tr)))
  s0_raw  <- run_ss(T_te[common_g, , drop = FALSE],  sig_list[pruned_sigs])
  s30_raw <- run_ss(mix30[common_g, , drop = FALSE], sig_list[pruned_sigs])
  ext_obj <- purify_ext(mix30, T_tr, N_tr, nhas_genes)
  s30_ext <- run_ss(ext_obj$pur_mat, sig_list[pruned_sigs])
  
  eval_m <- function(s0, sa, lbl) {
    sd0 <- apply(s0, 1, sd); med0 <- apply(s0, 1, median); meda <- apply(sa, 1, median)
    tibble(
      Direction = paste0(train_name, " -> ", test_name),
      Method    = lbl,
      Mean_Alpha_Hat = if (lbl == "External_S2") round(mean(ext_obj$alpha_hat), 3) else NA_real_,
      Abs_Drift = round(mean(abs(rowMeans(sa - s0) / sd0)), 3),
      Orig_Rho  = round(mean(row_spearman_fast(s0, sa)), 3),
      Fixed_Flip = round(mean(rowMeans((s0 >= med0) != (sa >= med0)) * 100), 1),
      Rel_Flip   = round(mean(rowMeans((s0 >= med0) != (sa >= meda)) * 100), 1)
    )
  }
  bind_rows(eval_m(s0_raw, s30_raw, "Raw_Unadjusted"),
            eval_m(s0_raw, s30_ext, "External_S2"))
}

ext_s2_df <- bind_rows(
  run_ext_s2_pair("ICGC_LIRI_JP", "TCGA_LIHC"),
  run_ext_s2_pair("TCGA_LIHC",    "ICGC_LIRI_JP")
)
write.csv(ext_s2_df, file.path(out_dir, "07a_external_S2_validation.csv"), row.names = FALSE)

message("\n>>> [结果 3] 严格跨队列外部参考 S2 验证 (Uniform 30% 混入):")
print(as.data.frame(ext_s2_df), row.names = FALSE)
message("\n>>> 07a 全部完成！结果已保存至 ", out_dir)
