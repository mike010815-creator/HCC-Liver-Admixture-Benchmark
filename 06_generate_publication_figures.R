# ==============================================================================
# 脚本名称: 06_generate_publication_figures.R
# 核心任务: 一键生成全部 6 张投稿级统一精修主图 (Figure 1 - Figure 6)
# 输出目录: 04_figures/final_manuscript/ (矢量 PDF + 300 DPI PNG + 300 DPI TIFF)
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Matrix)
  library(GSVA)
  library(Seurat)
  library(patchwork)
})

out_dir <- "04_figures/final_manuscript"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. 全局出版级主题、规范英文标签与统一色板
# ------------------------------------------------------------------------------
pretty_sig <- c(
  "CYT_Rooney"             = "Cytolytic Activity (CYT)",
  "T_cell_inflamed_GEP"    = "T-cell Inflamed GEP",
  "CD8_T_effector"         = "CD8+ T Effector",
  "HM_IFN_gamma"           = "Hallmark IFN-gamma",
  "HM_IFN_alpha"           = "Hallmark IFN-alpha",
  "HM_Allograft_Rejection" = "Hallmark Allograft Rejection",
  "Immune_Checkpoint"      = "Immune Checkpoint",
  "T_cell_exhaustion"      = "T-cell Exhaustion",
  "Treg_Signature"         = "Treg (CD4_c1)",
  "MDSC_Peng"              = "MDSC",
  "TAM_Peng"               = "TAM (Macrophage)",
  "HM_IL2_STAT5"           = "Hallmark IL2-STAT5",
  "HM_Complement"          = "Hallmark Complement",
  "HM_Inflammatory_Resp"   = "Hallmark Inflammatory Resp.",
  "HM_IL6_JAK_STAT3"       = "Hallmark IL6-JAK-STAT3",
  "HM_TNFA_via_NFKB"       = "Hallmark TNF-alpha/NF-kB",
  "Pan_F_TBRs"             = "Pan-F-TBRs (TGF-beta)",
  "CAF_Peng"               = "CAF",
  "WNT_target"             = "WNT Target",
  "HM_TGF_beta"            = "Hallmark TGF-beta",
  "HM_EMT"                 = "Hallmark EMT",
  "HM_Wnt_Beta_Catenin"    = "Hallmark Wnt/Beta-Catenin",
  "HM_Angiogenesis"        = "Hallmark Angiogenesis",
  "CTRL_G2M_Checkpoint"    = "G2M Checkpoint (Prolif. Ctrl)"
)

pretty_axis <- c(
  "1_Adaptive_Immunity"        = "Adaptive Immunity",
  "2_Suppression_Myeloid"      = "Suppression & Myeloid",
  "3_Innate_Inflammation"      = "Innate Inflammation",
  "4_Stromal_Exclusion"        = "Stromal & Exclusion",
  "0_Tumor_Proliferation_Ctrl" = "Proliferation Control"
)

axis_cols <- c(
  "Adaptive Immunity"     = "#3C5488",
  "Suppression & Myeloid" = "#00A087",
  "Innate Inflammation"   = "#E64B35",
  "Stromal & Exclusion"   = "#7E6148",
  "Proliferation Control" = "#555555"
)

pat_cols <- c(
  "Spurious Inflation"    = "#D62728",
  "Relatively Robust"     = "#2CA02C",
  "Proliferation Control" = "#2166AC"
)

map_sig  <- function(x) ifelse(x %in% names(pretty_sig),  unname(pretty_sig[x]),  x)
map_axis <- function(x) ifelse(x %in% names(pretty_axis), unname(pretty_axis[x]), x)

theme_pub <- function(base_size = 9.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.35),
      panel.border      = element_rect(color = "grey25", linewidth = 0.55),
      strip.background  = element_rect(fill = "#F3F5F8", color = "grey30", linewidth = 0.45),
      strip.text        = element_text(face = "bold", size = base_size - 0.5, color = "black"),
      axis.text         = element_text(color = "black", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size, face = "bold"),
      plot.title        = element_text(face = "bold", size = base_size + 1, hjust = 0, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = base_size - 1.2, color = "grey35", margin = margin(b = 5)),
      legend.title      = element_text(face = "bold", size = base_size - 0.8),
      legend.text       = element_text(size = base_size - 1.2),
      legend.key.size   = unit(0.42, "cm"),
      plot.margin       = margin(6, 10, 6, 8)
    )
}

save_pub_fig <- function(plot_obj, stem, width, height) {
  p_pdf  <- file.path(out_dir, paste0(stem, ".pdf"))
  p_png  <- file.path(out_dir, paste0(stem, ".png"))
  p_tiff <- file.path(out_dir, paste0(stem, ".tiff"))
  
  ggsave(p_pdf, plot_obj, width = width, height = height, device = "pdf")
  ggsave(p_png, plot_obj, width = width, height = height, dpi = 300, bg = "white")
  tryCatch(
    ggsave(p_tiff, plot_obj, width = width, height = height, dpi = 300,
           device = "tiff", compression = "lzw", bg = "white"),
    error = function(e) message("    (TIFF 跳过: ", e$message, ")")
  )
  message(sprintf("    已导出 -> %s (.pdf / .png / .tiff)", stem))
}

# ------------------------------------------------------------------------------
# 加载全项目已生成的清洗对象
# ------------------------------------------------------------------------------
tcga_obj <- readRDS("01_clean_data/tcga_lihc_clean.rds")
toil_obj <- readRDS("01_clean_data/toil_lihc_gtex_clean.rds")
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
sim_all  <- readRDS("01_clean_data/02d_multicohort_sim_results.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
res3_cal <- readRDS("01_clean_data/03b_scrna_calibrated.rds")
res4_cli <- readRDS("01_clean_data/04b_clinical_survival_results.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

sig_list   <- sigs_obj$sig_list
sig_meta   <- sigs_obj$sig_meta
nhas_genes <- res3_obj$nhas_genes

# ==============================================================================
# FIGURE 1: TCGA-LIHC 24 特征异质性敏感性全景与三大打分范式对比
# ==============================================================================
message(">>> [1/6] 正在生成投稿级 Figure 1（统一 24 特征并补齐 GSVA 对比）...")

# 补算 24 特征在 TCGA 30% 混入下的 GSVA Spearman rho (耗时 ~1 秒)
m0_tcga  <- log2(tcga_obj$T_tpm + 1)
m30_tcga <- log2(0.7 * tcga_obj$T_tpm + 0.3 * tcga_obj$N_tpm + 1)
g0  <- gsva(gsvaParam(as.matrix(m0_tcga),  sig_list, minSize = 2, kcdf = "Gaussian"), verbose = FALSE)
g30 <- gsva(gsvaParam(as.matrix(m30_tcga), sig_list, minSize = 2, kcdf = "Gaussian"), verbose = FALSE)

gsva_rho30 <- tibble(
  Signature    = rownames(g0),
  Method       = "GSVA",
  Spearman_Rho = sapply(rownames(g0), function(s) cor(as.numeric(g0[s, ]), as.numeric(g30[s, ]), method = "spearman"))
)

tcga_ss <- sim_all %>%
  filter(grepl("TCGA", Cohort), Method == "ssGSEA") %>%
  mutate(Sig_Label = map_sig(Signature), Axis_Label = map_axis(Axis))

tcga_sum30 <- tcga_ss %>%
  filter(Alpha_Pct == 30) %>%
  mutate(Pattern = case_when(
    Signature == "CTRL_G2M_Checkpoint" ~ "Proliferation Control",
    Std_Drift > 0.30                   ~ "Spurious Inflation",
    TRUE                               ~ "Relatively Robust"
  )) %>%
  arrange(desc(Std_Drift))

sig_ord1 <- rev(tcga_sum30$Sig_Label)

# Panel A: 4 大免疫轴 + 增殖对照的 0%–50% 轨迹曲线 (1x4 横排更紧凑省空间)
p1a_df <- tcga_ss %>%
  left_join(tcga_sum30 %>% select(Signature, Pattern), by = "Signature") %>%
  mutate(Facet_Axis = ifelse(Axis_Label == "Proliferation Control", "Adaptive Immunity", Axis_Label),
         Facet_Axis = factor(Facet_Axis, levels = c("Adaptive Immunity", "Suppression & Myeloid",
                                                    "Innate Inflammation", "Stromal & Exclusion")))

p1a <- ggplot(p1a_df, aes(x = Alpha_Pct, y = Std_Drift, group = Signature,
                          color = Pattern, linetype = Pattern)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_hline(yintercept = 0.3, linetype = "dotted", color = "grey65", linewidth = 0.45) +
  geom_line(linewidth = 0.85, alpha = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~ Facet_Axis, nrow = 1) +
  scale_color_manual(values = pat_cols, name = "Drift Pattern (at 30% Admixture)") +
  scale_linetype_manual(values = c("Spurious Inflation" = "solid",
                                   "Relatively Robust"  = "solid",
                                   "Proliferation Control" = "dashed"),
                        name = "Drift Pattern (at 30% Admixture)") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50)) +
  labs(title = "Score Drift Trajectories Across Non-tumor Admixture Gradients (TCGA-LIHC, n=50 pairs)",
       x = "Non-tumor Liver Admixture (%)", y = "Standardized Drift (SD)") +
  theme_pub(9.5) +
  theme(legend.position = "bottom")

# Panel B: 30% 混入标准化偏移棒棒糖图
p1b <- tcga_sum30 %>%
  mutate(Sig_Label = factor(Sig_Label, levels = sig_ord1)) %>%
  ggplot(aes(x = Std_Drift, y = Sig_Label, color = Pattern)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = Std_Drift, yend = Sig_Label), linewidth = 1.0) +
  geom_point(size = 2.6) +
  scale_color_manual(values = pat_cols, guide = "none") +
  labs(title = "Score Drift at 30%", x = "Standardized Drift (SD)", y = NULL) +
  theme_pub(9.5)

# Panel C: 5%–50% 固定阈值临床分层翻转率热图
p1c <- tcga_ss %>%
  filter(Alpha_Pct > 0) %>%
  mutate(Sig_Label = factor(Sig_Label, levels = sig_ord1),
         Alpha_Fac = factor(paste0(Alpha_Pct, "%"),
                            levels = paste0(c(5, 10, 20, 30, 40, 50), "%"))) %>%
  ggplot(aes(x = Alpha_Fac, y = Sig_Label, fill = Flip_Fixed_Pct)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f", Flip_Fixed_Pct)), size = 2.6, color = "black") +
  scale_fill_gradient2(low = "#F7FBFF", mid = "#FDAE6B", high = "#BD0026",
                       midpoint = 22, name = "Flip Rate (%)") +
  labs(title = "Misclassification Rate (%)", x = "Admixture Proportion", y = NULL) +
  theme_pub(9.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.key.width = unit(0.9, "cm"))

# Panel D: 三大打分算法在 30% 混入下的秩次保真度对比
rho30_all <- bind_rows(
  sim_all %>% filter(grepl("TCGA", Cohort), Alpha_Pct == 30) %>% select(Signature, Method, Spearman_Rho),
  gsva_rho30
) %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = sig_ord1),
    Method    = factor(Method, levels = c("GSVA", "ssGSEA", "Zscore"),
                       labels = c("GSVA", "ssGSEA", "Z-score"))
  )

p1d <- ggplot(rho30_all, aes(x = Method, y = Sig_Label, fill = Spearman_Rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Spearman_Rho)), size = 2.6, color = "black") +
  scale_fill_gradient2(low = "#D73027", mid = "#FEE08B", high = "#1A9850",
                       midpoint = 0.85, limits = c(0.60, 1.00), oob = scales::squish,
                       name = "Spearman rho") +
  labs(title = "Rank Stability (30%)", x = "Scoring Paradigm", y = NULL) +
  theme_pub(9.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.key.width = unit(0.9, "cm"))

fig1 <- (p1a / (p1b | p1c | p1d + plot_layout(widths = c(1.55, 1.35, 0.95)))) +
  plot_layout(heights = c(1.0, 2.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig1, "Figure1_Admixture_Sensitivity_Landscape", width = 14.5, height = 10.5)

# ==============================================================================
# FIGURE 2: 混入来源异质性 (癌旁肝 vs GTEx 健康肝) 与三大独立队列跨平台复现
# ==============================================================================
message(">>> [2/6] 正在生成投稿级 Figure 2（重算 24 特征癌旁 vs GTEx 并修复图例排版）...")

# 用更新后的完整 24 特征库快速重算 Toil 队列 30% 癌旁 vs GTEx 健康肝 (ssGSEA)
run_ss_toil <- function(lin_mat) {
  m <- log2(lin_mat + 1)
  gsva(ssgseaParam(as.matrix(m), sig_list, minSize = 2, normalize = FALSE), verbose = FALSE)
}
s0_toil  <- run_ss_toil(toil_obj$T_toil)
s_adj30  <- run_ss_toil(0.7 * toil_obj$T_toil + 0.3 * toil_obj$N_adj)

drift_adj <- tibble(
  Signature = rownames(s0_toil),
  Source    = "TCGA Adjacent Liver (30%)",
  Std_Drift = sapply(rownames(s0_toil), function(s) mean(s_adj30[s, ] - s0_toil[s, ]) / sd(s0_toil[s, ]))
)

n_toil <- ncol(toil_obj$T_toil)
drift_gtex <- map_dfr(1:5, function(rep_i) {
  set.seed(2026 + rep_i)
  pick <- sample(ncol(toil_obj$N_gtex), n_toil, replace = FALSE)
  sm   <- run_ss_toil(0.7 * toil_obj$T_toil + 0.3 * toil_obj$N_gtex[, pick, drop = FALSE])
  tibble(
    Signature = rownames(s0_toil),
    Std_Drift = sapply(rownames(s0_toil), function(s) mean(sm[s, ] - s0_toil[s, ]) / sd(s0_toil[s, ]))
  )
}) %>%
  group_by(Signature) %>%
  summarise(Std_Drift = mean(Std_Drift), .groups = "drop") %>%
  mutate(Source = "GTEx Healthy Liver (30%)")

toil_cmp <- bind_rows(drift_adj, drift_gtex) %>%
  mutate(Sig_Label = map_sig(Signature))

# 统一排序：按三队列 ssGSEA 30% TCGA 漂移从高到低
wide30 <- sim_all %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  mutate(C_short = case_when(grepl("TCGA", Cohort) ~ "TCGA",
                             grepl("ICGC", Cohort) ~ "ICGC",
                             TRUE                  ~ "GSE14520")) %>%
  select(Signature, Axis, C_short, Std_Drift) %>%
  pivot_wider(names_from = C_short, values_from = Std_Drift) %>%
  mutate(Sig_Label = map_sig(Signature), Axis_Label = map_axis(Axis)) %>%
  arrange(desc(TCGA))

sig_ord2 <- rev(wide30$Sig_Label)

# Panel A: 癌旁肝 vs GTEx 健康肝哑铃对比图 (24 特征全齐)
p2a <- toil_cmp %>%
  mutate(Sig_Label = factor(Sig_Label, levels = sig_ord2)) %>%
  ggplot(aes(x = Std_Drift, y = Sig_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey65", linewidth = 0.85) +
  geom_point(aes(color = Source, shape = Source), size = 2.8) +
  scale_color_manual(values = c("TCGA Adjacent Liver (30%)" = "#D62728",
                                "GTEx Healthy Liver (30%)"  = "#1F77B4"),
                     name = NULL) +
  scale_shape_manual(values = c("TCGA Adjacent Liver (30%)" = 16,
                                "GTEx Healthy Liver (30%)"  = 17),
                     name = NULL) +
  labs(title = "Adjacent Inflamed vs. Healthy Liver",
       subtitle = "Toil-harmonized cohort (ssGSEA, 30%)",
       x = "Standardized Drift (SD)", y = NULL) +
  theme_pub(9.5) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# Panel B: 三大独立队列 30% 混入热图
p2b <- sim_all %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = sig_ord2),
    C_lbl = case_when(
      grepl("TCGA", Cohort) ~ "TCGA-LIHC\n(n=50, RNA-seq)",
      grepl("ICGC", Cohort) ~ "ICGC-LIRI-JP\n(n=175, RNA-seq)",
      TRUE                  ~ "GSE14520\n(n=214, Array)"
    ),
    C_lbl = factor(C_lbl, levels = c("TCGA-LIHC\n(n=50, RNA-seq)",
                                     "ICGC-LIRI-JP\n(n=175, RNA-seq)",
                                     "GSE14520\n(n=214, Array)"))
  ) %>%
  ggplot(aes(x = C_lbl, y = Sig_Label, fill = Std_Drift)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", Std_Drift)), size = 2.7, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Std Drift (30%)") +
  labs(title = "Cross-Cohort Reproducibility",
       subtitle = "439 paired HCC-adjacent samples",
       x = NULL, y = NULL) +
  theme_pub(9.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 8.2),
        legend.position = "bottom", legend.key.width = unit(1.0, "cm"))

# Panel C & D: 跨队列散点图 (彻底修复图例换行不裁切)
r_ti <- cor(wide30$TCGA, wide30$ICGC, method = "spearman")
r_tg <- cor(wide30$TCGA, wide30$GSE14520, method = "spearman")

p2c <- ggplot(wide30, aes(x = TCGA, y = ICGC, color = Axis_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7, fill = "grey88") +
  geom_point(size = 2.6, alpha = 0.92) +
  scale_color_manual(values = axis_cols, name = "Biological Axis") +
  annotate("text", x = min(wide30$TCGA), y = max(wide30$ICGC), hjust = 0, vjust = 1,
           fontface = "bold", size = 3.2,
           label = sprintf("Spearman rho = %.3f\n(RNA-seq vs. RNA-seq)", r_ti)) +
  labs(title = "TCGA-LIHC vs. ICGC-LIRI-JP",
       x = "TCGA-LIHC 30% Drift (SD)", y = "ICGC-LIRI-JP 30% Drift (SD)") +
  theme_pub(9.5) +
  theme(legend.position = "none")

p2d <- ggplot(wide30, aes(x = TCGA, y = GSE14520, color = Axis_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.7, fill = "grey88") +
  geom_point(size = 2.6, alpha = 0.92) +
  scale_color_manual(values = axis_cols, name = "Biological Axis") +
  annotate("text", x = min(wide30$TCGA), y = max(wide30$GSE14520), hjust = 0, vjust = 1,
           fontface = "bold", size = 3.2,
           label = sprintf("Spearman rho = %.3f\n(RNA-seq vs. Microarray)", r_tg)) +
  labs(title = "TCGA-LIHC vs. GSE14520",
       x = "TCGA-LIHC 30% Drift (SD)", y = "GSE14520 30% Drift (SD)") +
  theme_pub(9.5) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = "bottom")

fig2 <- (p2a | p2b | (p2c / p2d)) +
  plot_layout(widths = c(1.25, 1.05, 1.10)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig2, "Figure2_Source_and_MultiCohort_Validation", width = 15, height = 8.6)

# ==============================================================================
# FIGURE 3: 单细胞机制拆解 (生理比例校正 + 峰值验证 + 癌旁细胞来源溯源)
# ==============================================================================
message(">>> [3/6] 正在生成投稿级 Figure 3...")

tri_wide <- res3_cal$tri_wide %>% mutate(Sig_Label = map_sig(Signature))
sig_ord3 <- rev(tri_wide$Sig_Label)

p3a <- res3_cal$sweep_df %>%
  filter(Hep_Frac %in% c(100, 80, 0)) %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = sig_ord3),
    Arm_Label = factor(
      case_when(
        Hep_Frac == 100 ~ "Pure Hepatocytes (100% Hep)",
        Hep_Frac == 80  ~ "Physiological Liver (80% Hep + 20% NonHep)",
        Hep_Frac == 0   ~ "Pure Non-Parenchymal (0% Hep)"
      ),
      levels = c("Pure Hepatocytes (100% Hep)",
                 "Physiological Liver (80% Hep + 20% NonHep)",
                 "Pure Non-Parenchymal (0% Hep)")
    )
  ) %>%
  ggplot(aes(x = Std_Drift, y = Sig_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey72", linewidth = 0.85) +
  geom_point(aes(color = Arm_Label, shape = Arm_Label), size = 2.8) +
  scale_color_manual(values = c("#1F77B4", "#D62728", "#FF7F0E"), name = "Admixture Composition") +
  scale_shape_manual(values = c(16, 17, 15), name = "Admixture Composition") +
  labs(title = "Calibrated Single-Cell Mechanism Decomposition",
       subtitle = "Paired scRNA-seq pseudo-bulk (GSE149614, n=7 patients, 30% admixture)",
       x = "Standardized Drift (ssGSEA, 30%)", y = NULL) +
  theme_pub(9.5) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p3b <- res3_cal$sweep_cor %>%
  pivot_longer(cols = c(Rho_vs_TCGA, Rho_vs_ICGC),
               names_to = "Bulk_Cohort", values_to = "Spearman_Rho") %>%
  mutate(Bulk_Cohort = recode(Bulk_Cohort,
                              Rho_vs_TCGA = "vs. TCGA-LIHC Bulk",
                              Rho_vs_ICGC = "vs. ICGC-LIRI-JP Bulk")) %>%
  ggplot(aes(x = Hep_Frac, y = Spearman_Rho, color = Bulk_Cohort)) +
  annotate("rect", xmin = 70, xmax = 85, ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.7) +
  annotate("text", x = 77.5, y = 0.24, label = "Physiological\nliver mRNA\n(70-85%)",
           size = 2.9, color = "#1B5E20", fontface = "bold") +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("vs. TCGA-LIHC Bulk" = "#B2182B",
                                "vs. ICGC-LIRI-JP Bulk" = "#2166AC"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  labs(title = "Correction of scRNA-seq Dissociation Bias",
       subtitle = "Peak bulk concordance at physiological hepatocyte fraction",
       x = "Hepatocyte RNA Fraction in Adjacent Liver (%)",
       y = "Concordance with Bulk Drift (Spearman rho)") +
  theme_pub(9.5) +
  theme(legend.position = "bottom")

# Panel C: 污染易感基因的癌旁主导细胞来源堆叠条形图
cell_cols <- c("T/NK" = "#E64B35", "Myeloid" = "#00A087", "Endothelial" = "#4DBBD5",
               "Fibroblast" = "#3C5488", "Hepatocyte" = "#F39B7F", "B" = "#8491B4")

p3c <- res3_obj$sig_gene_map %>%
  filter(Admixture_Prone) %>%
  mutate(
    Sig_Label  = factor(map_sig(Signature), levels = sig_ord3),
    Cell_Group = sub("^N_", "", Top_NTL_CellType)
  ) %>%
  count(Sig_Label, Cell_Group) %>%
  ggplot(aes(x = n, y = Sig_Label, fill = Cell_Group)) +
  geom_col(width = 0.72) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = cell_cols, name = "Dominant NTL Cell Source") +
  labs(title = "Cellular Origin of Admixture-Prone Genes",
       subtitle = "Adjacent-enriched genes (log2FC > 0.8, CPM > 10)",
       x = "Number of Admixture-Prone Genes", y = NULL) +
  theme_pub(9.5) +
  guides(fill = guide_legend(nrow = 2)) +
  theme(axis.text.y = element_text(size = 7.8), legend.position = "bottom")

fig3 <- (p3a | (p3b / p3c + plot_layout(heights = c(0.85, 1.35)))) +
  plot_layout(widths = c(1.15, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig3, "Figure3_scRNA_Mechanism_Decomposition", width = 14.5, height = 9.5)

# ==============================================================================
# FIGURE 4: 去 UMI 混杂空间转录组验证 (HCC3R 切片分布 + 跨模态验证 + 空间净化)
# ==============================================================================
message(">>> [4/6] 正在生成投稿级 Figure 4（读取 HCC3R 切片与去混杂结果）...")

sp_clean_df <- fread("03_results/04c_hcc3r_spatial_deconfounded.csv", data.table = FALSE) %>%
  mutate(Sig_Label = map_sig(Signature))
h3_dir <- "01_clean_data/spatial_unpacked/GSE238264/HCC3R"
stopifnot(dir.exists(h3_dir))

cnt3 <- Read10X(file.path(h3_dir, "filtered_feature_bc_matrix"))
if (is.list(cnt3)) cnt3 <- cnt3[[1]]
seu3 <- CreateSeuratObject(counts = cnt3, assay = "Spatial", project = "HCC3R")
img3 <- Read10X_Image(file.path(h3_dir, "spatial"), image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
DefaultAssay(img3) <- "Spatial"
common_bc <- intersect(colnames(seu3), Cells(img3))
seu3 <- subset(seu3, cells = common_bc)
img3 <- img3[common_bc]
seu3[["slice1"]] <- img3
seu3 <- subset(seu3, subset = nFeature_Spatial >= 500)
seu3 <- NormalizeData(seu3, verbose = FALSE)

score_ms <- function(obj, gvec, tag) {
  gu <- intersect(gvec, rownames(obj))
  tmp <- AddModuleScore(obj, features = list(gu), name = tag, ctrl = 50, seed = 2026)
  as.numeric(tmp[[paste0(tag, "1")]][, 1])
}

seu3$SP_NHAS     <- score_ms(seu3, nhas_genes, "M_NHAS_")
seu3$SP_G2M      <- score_ms(seu3, sig_list$CTRL_G2M_Checkpoint, "M_G2M_")
seu3$SP_Comp     <- score_ms(seu3, sig_list$HM_Complement, "M_COMP_")
seu3$SP_Comp_Pur <- score_ms(seu3, res5_obj$purified_sig_list$HM_Complement, "M_COMPP_")
seu3$SP_TNFA     <- score_ms(seu3, sig_list$HM_TNFA_via_NFKB, "M_TNFA_")
seu3$SP_TNFA_Pur <- score_ms(seu3, res5_obj$purified_sig_list$HM_TNFA_via_NFKB, "M_TNFAP_")

# 展示最具病理对比意义的 6 张空间分布图 (镜像对照 + 补体与 TNFa 净化前后对比)
sp_specs <- c(
  "SP_NHAS"     = "NHAS (Non-tumor Hepatocyte)",
  "SP_G2M"      = "G2M Checkpoint (Prolif. Ctrl)",
  "SP_Comp"     = "Hallmark Complement (Raw)",
  "SP_Comp_Pur" = "Hallmark Complement (Purified)",
  "SP_TNFA"     = "Hallmark TNF-alpha/NF-kB (Raw)",
  "SP_TNFA_Pur" = "Hallmark TNF-alpha/NF-kB (Purified)"
)

sp_maps <- lapply(names(sp_specs), function(cn) {
  med_v <- median(seu3[[cn]][, 1], na.rm = TRUE)
  suppressMessages(
    SpatialFeaturePlot(seu3, features = cn, image.alpha = 0.35, pt.size.factor = 1.65) +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = med_v, name = "Score") +
      ggtitle(sp_specs[[cn]]) +
      theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
            legend.position = "right",
            legend.key.height = unit(0.4, "cm"),
            legend.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 7.5),
            plot.margin = margin(2, 4, 2, 4))
  )
})
p4a <- wrap_plots(sp_maps, ncol = 3)

# Panel B: 24 特征去混杂空间共变条形图
r_sp_bulk <- cor(sp_clean_df$Spatial_MS_Rho_Raw, sp_clean_df$Bulk30_Mean, method = "spearman")
sig_ord4  <- rev(sp_clean_df$Sig_Label)

p4b <- sp_clean_df %>%
  mutate(
    Sig_Label = factor(Sig_Label, levels = sig_ord4),
    Dir = ifelse(Spatial_MS_Rho_Raw >= 0, "Co-enriched with NHAS (Positive)", "Tumor-Core Enriched (Negative)")
  ) %>%
  ggplot(aes(x = Spatial_MS_Rho_Raw, y = Sig_Label, fill = Dir)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
  geom_col(width = 0.72, alpha = 0.9) +
  scale_fill_manual(values = c("Co-enriched with NHAS (Positive)" = "#D62728",
                               "Tumor-Core Enriched (Negative)"   = "#2166AC"), name = NULL) +
  labs(title = "Deconfounded Spatial Co-variation (HCC3R)",
       subtitle = sprintf("Concordance with 3-cohort Bulk 30%% drift: Spearman rho = %.3f", r_sp_bulk),
       x = "Spatial Spearman rho with NHAS", y = NULL) +
  theme_pub(9.5) +
  theme(legend.position = "bottom")

# Panel C: 空间切片上 S3 基因级净化的去共变效果 (突出展示正向共变最强的炎症/补体通路)
sp_pur_df <- sp_clean_df %>% filter(!is.na(Spatial_MS_Rho_Pur))
pur_ord4  <- sp_pur_df %>% arrange(Spatial_MS_Rho_Raw) %>% pull(Sig_Label)
w_sp      <- wilcox.test(sp_pur_df$Spatial_MS_Rho_Raw, sp_pur_df$Spatial_MS_Rho_Pur, paired = TRUE)

p4c <- sp_pur_df %>%
  select(Sig_Label, Spatial_MS_Rho_Raw, Spatial_MS_Rho_Pur) %>%
  pivot_longer(-Sig_Label, names_to = "Ver", values_to = "Rho") %>%
  mutate(
    Sig_Label = factor(Sig_Label, levels = pur_ord4),
    Ver = recode(Ver, Spatial_MS_Rho_Raw = "Original Signature (Raw)",
                 Spatial_MS_Rho_Pur = "scRNA-Purified Signature (S3)")
  ) %>%
  ggplot(aes(x = Rho, y = Sig_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey65", linewidth = 0.85) +
  geom_point(aes(color = Ver, shape = Ver), size = 2.8) +
  scale_color_manual(values = c("Original Signature (Raw)"      = "#D62728",
                                "scRNA-Purified Signature (S3)" = "#00A087"), name = NULL) +
  scale_shape_manual(values = c("Original Signature (Raw)" = 16,
                                "scRNA-Purified Signature (S3)" = 17), name = NULL) +
  labs(title = "Spatial Deconfounding via scRNA Gene Purification",
       subtitle = sprintf("18 pruned signatures (paired Wilcoxon P = %.3f; strongest drop in innate inflammation)",
                          w_sp$p.value),
       x = "Spatial Spearman rho with NHAS (HCC3R)", y = NULL) +
  theme_pub(9.5) +
  theme(legend.position = "bottom")

fig4 <- (p4a / (p4b | p4c)) +
  plot_layout(heights = c(1.05, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig4, "Figure4_Spatial_Validation", width = 14.5, height = 13.0)

# ==============================================================================
# FIGURE 5: 789 例真实初治原发 HCC 的共变验证与临床预后偏倚 (新增分层森林图 Panel D)
# ==============================================================================
message(">>> [5/6] 正在生成投稿级 Figure 5（含 Low vs High-NHAS 分层森林图）...")

real_cor_wide <- res4_cli$real_cor_wide %>%
  mutate(Sig_Label = map_sig(Signature), Axis_Label = map_axis(Axis))
sig_ord5 <- rev(real_cor_wide$Sig_Label)

# Panel A: 三大真实队列未混合原发癌 NHAS 相关热图
p5a <- res4_cli$real_cor_df %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = sig_ord5),
    Cohort_lbl = factor(Cohort,
                        levels = c("TCGA_LIHC", "ICGC_LIRI_JP", "GSE14520"),
                        labels = c("TCGA-LIHC\n(n=365)", "ICGC-LIRI-JP\n(n=203)", "GSE14520\n(n=221)")),
    star = ifelse(Real_P < 0.001, "***", ifelse(Real_P < 0.01, "**", ifelse(Real_P < 0.05, "*", "")))
  ) %>%
  ggplot(aes(x = Cohort_lbl, y = Sig_Label, fill = Real_Rho_with_NHAS)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f%s", Real_Rho_with_NHAS, star)), size = 2.65) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Spearman rho\nwith NHAS") +
  labs(title = "Real-Tumor Co-variation with NHAS",
       subtitle = "789 unmixed treatment-naive primary HCCs", x = NULL, y = NULL) +
  theme_pub(9.5) +
  theme(axis.text.x = element_text(face = "bold", size = 8.2),
        legend.position = "bottom", legend.key.width = unit(1.0, "cm"))

# Panel B: 模拟预测 vs 真实观察散点图
c_mean <- cor(real_cor_wide$Mean_Real_Rho, real_cor_wide$Mean_Sim30_Drift, method = "spearman")
c_icgc <- cor(real_cor_wide$RealRho_ICGC_LIRI_JP, real_cor_wide$Sim30_ICGC, method = "spearman")

p5b <- ggplot(real_cor_wide, aes(x = Mean_Sim30_Drift, y = Mean_Real_Rho, color = Axis_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey25", linewidth = 0.75, fill = "grey88") +
  geom_point(size = 2.7, alpha = 0.92) +
  scale_color_manual(values = axis_cols, name = "Biological Axis") +
  annotate("text", x = min(real_cor_wide$Mean_Sim30_Drift), y = max(real_cor_wide$Mean_Real_Rho),
           hjust = 0, vjust = 1, fontface = "bold", size = 3.1,
           label = sprintf("Pooled 3-Cohort: rho = %.3f\nICGC-LIRI-JP: rho = %.3f", c_mean, c_icgc)) +
  labs(title = "In Silico vs. Real-World Observation",
       subtitle = "Admixture-differentiation coupling in real tumors",
       x = "Predicted 30% Admixture Drift (Mean SD)",
       y = "Observed Real-Tumor NHAS Correlation (Mean rho)") +
  theme_pub(9.5) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = "bottom")

# Panel C: NHAS 校正前 vs 校正后 Meta-Cox HR 偏移
meta_cox <- res4_cli$meta_cox_df %>% mutate(Sig_Label = map_sig(Signature))
p5c <- meta_cox %>%
  select(Sig_Label, HR_Raw, HR_Adj) %>%
  pivot_longer(c(HR_Raw, HR_Adj), names_to = "Model", values_to = "HR") %>%
  mutate(
    Sig_Label = factor(Sig_Label, levels = sig_ord5),
    Model = recode(Model, HR_Raw = "Unadjusted Cox (~ Sig)",
                   HR_Adj = "NHAS-Adjusted Cox (~ Sig + NHAS)")
  ) %>%
  ggplot(aes(x = HR, y = Sig_Label)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey65", linewidth = 0.85) +
  geom_point(aes(color = Model, shape = Model), size = 2.7) +
  scale_color_manual(values = c("Unadjusted Cox (~ Sig)" = "#D62728",
                                "NHAS-Adjusted Cox (~ Sig + NHAS)" = "#1F77B4"), name = NULL) +
  scale_shape_manual(values = c(17, 16), name = NULL) +
  labs(title = "NHAS Confounding in Meta-Cox (n=789)",
       subtitle = "Pooled OS HR per 1-SD increase",
       x = "Pooled Hazard Ratio (OS)", y = NULL) +
  theme_pub(9.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "vertical")

# Panel D (重磅新增): Low-NHAS vs High-NHAS 分层 Meta-Cox 森林图 (精选 10 个临床关键/显著交互特征)
key_strat <- c("CYT_Rooney", "CD8_T_effector", "T_cell_inflamed_GEP", "HM_IFN_gamma", "HM_IFN_alpha",
               "Pan_F_TBRs", "CAF_Peng", "HM_EMT", "MDSC_Peng", "CTRL_G2M_Checkpoint")
p5d_df <- meta_cox %>%
  filter(Signature %in% key_strat) %>%
  mutate(
    Sig_Star = ifelse(P_Int_Meta < 0.01, paste0(Sig_Label, " **"),
                      ifelse(P_Int_Meta < 0.05, paste0(Sig_Label, " *"), Sig_Label))
  )

ord_5d <- rev(p5d_df$Sig_Star[match(key_strat, p5d_df$Signature)])

p5d_long <- bind_rows(
  p5d_df %>% transmute(Sig_Star, Stratum = "Low-NHAS (Purer Tumors)", HR = HR_Low, Lo = Lo_Low, Hi = Hi_Low, P = P_Low),
  p5d_df %>% transmute(Sig_Star, Stratum = "High-NHAS (Higher Admixture)", HR = HR_High, Lo = Lo_High, Hi = Hi_High, P = P_High)
) %>%
  mutate(
    Sig_Star = factor(Sig_Star, levels = ord_5d),
    Stratum  = factor(Stratum, levels = c("Low-NHAS (Purer Tumors)", "High-NHAS (Higher Admixture)")),
    Sig_Flag = ifelse(P < 0.05, "P < 0.05", "NS")
  )

p5d <- ggplot(p5d_long, aes(x = HR, y = Sig_Star, color = Stratum)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = Lo, xmax = Hi), height = 0.22, linewidth = 0.7,
                 position = position_dodge(width = 0.52)) +
  geom_point(aes(shape = Sig_Flag), size = 2.6, position = position_dodge(width = 0.52)) +
  scale_color_manual(values = c("Low-NHAS (Purer Tumors)"      = "#1F77B4",
                                "High-NHAS (Higher Admixture)" = "#D62728"), name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "NS" = 1), name = NULL) +
  labs(title = "NHAS-Stratified Prognostic Attenuation",
       subtitle = "3-cohort Meta-HR (95% CI); * P_int < 0.05, ** P_int < 0.01",
       x = "Pooled Hazard Ratio (95% CI)", y = NULL) +
  theme_pub(9.5) +
  theme(legend.position = "bottom", legend.direction = "vertical")

fig5 <- ((p5a | p5c) | (p5b / p5d + plot_layout(heights = c(0.95, 1.15)))) +
  plot_layout(widths = c(1.35, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig5, "Figure5_Clinical_Prognostic_Bias", width = 15.5, height = 9.5)

# ==============================================================================
# FIGURE 6: 三大纠偏策略基准测试 + 单细胞净化降漂 + 789 例临床预后“解封”验证
# ==============================================================================
message(">>> [6/6] 正在生成投稿级 Figure 6（修复纵轴裁切并加入 Panel E 临床解封验证）...")

bench_df    <- res5_obj$bench_df
pruned_sigs <- res5_obj$prune_summary %>% filter(Was_Pruned) %>% pull(Signature)
fidelity_df <- res5_obj$fidelity_df

strat_labels <- c(
  "0_Unadjusted_Raw"          = "Unadjusted\n(Raw)",
  "1_NHAS_Residualization"    = "S1: NHAS\nResidual",
  "2_Linear_Ref_Subtraction"  = "S2: Linear\nPurification",
  "3_scRNA_Gene_Purification" = "S3: scRNA\nPurification"
)
strat_cols <- c(
  "0_Unadjusted_Raw"          = "#D62728",
  "1_NHAS_Residualization"    = "#FF7F0E",
  "2_Linear_Ref_Subtraction"  = "#1F77B4",
  "3_scRNA_Gene_Purification" = "#00A087"
)

# Panel A: 均匀 30% 混入绝对漂移 (纵轴标签简洁清晰，彻底杜绝裁切)
p6a <- bench_df %>%
  filter(Scenario == "Uniform_30pct", Signature %in% pruned_sigs) %>%
  ggplot(aes(x = Strategy, y = Abs_Drift, fill = Strategy)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.9, linewidth = 0.45) +
  scale_x_discrete(labels = strat_labels) +
  scale_fill_manual(values = strat_cols, guide = "none") +
  labs(title = "Absolute Drift (Uniform 30%)", x = NULL, y = "|Standardized Drift| (SD)") +
  theme_pub(9.2)

# Panel B: 异质性 0%–40% 混入秩次恢复
p6b <- bench_df %>%
  filter(Scenario == "Heterogeneous_0_40pct", Signature %in% pruned_sigs) %>%
  ggplot(aes(x = Strategy, y = Spearman_Rho, fill = Strategy)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.9, linewidth = 0.45) +
  scale_x_discrete(labels = strat_labels) +
  scale_fill_manual(values = strat_cols, guide = "none") +
  labs(title = "Rank Stability (Variable 0-40%)", x = NULL, y = "Spearman rho vs. Pure") +
  theme_pub(9.2)

# Panel C: 异质性 0%–40% 混入分层翻转率
p6c <- bench_df %>%
  filter(Scenario == "Heterogeneous_0_40pct", Signature %in% pruned_sigs) %>%
  ggplot(aes(x = Strategy, y = Flip_Fixed_Pct, fill = Strategy)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.9, linewidth = 0.45) +
  scale_x_discrete(labels = strat_labels) +
  scale_fill_manual(values = strat_cols, guide = "none") +
  labs(title = "Misclassification (Variable 0-40%)", x = NULL, y = "Flip Rate (%)") +
  theme_pub(9.2)

# Panel D: 18 个净化特征在三队列 30% 混入下的漂移压制 + 右侧标注纯肿瘤保真度 rho
s3_df <- bench_df %>%
  filter(Scenario == "Uniform_30pct",
         Strategy %in% c("0_Unadjusted_Raw", "3_scRNA_Gene_Purification"),
         Signature %in% pruned_sigs) %>%
  group_by(Signature, Strategy) %>%
  summarise(Mean_Drift = mean(Std_Drift), .groups = "drop")

s3_ord_df <- s3_df %>%
  pivot_wider(names_from = Strategy, values_from = Mean_Drift) %>%
  mutate(Red = abs(`0_Unadjusted_Raw`) - abs(`3_scRNA_Gene_Purification`),
         Sig_Label = map_sig(Signature)) %>%
  left_join(fidelity_df, by = "Signature") %>%
  arrange(desc(Red))

s3_ord6 <- rev(s3_ord_df$Sig_Label)

p6d <- s3_df %>%
  mutate(Sig_Label = factor(map_sig(Signature), levels = s3_ord6)) %>%
  ggplot(aes(x = Mean_Drift, y = Sig_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey65", linewidth = 0.9) +
  geom_point(aes(color = Strategy, shape = Strategy), size = 2.8) +
  geom_text(
    data = s3_ord_df %>% mutate(Sig_Label = factor(Sig_Label, levels = s3_ord6)),
    aes(x = 0.78, y = Sig_Label, label = sprintf("rho=%.2f", Mean_PureTumor_Fidelity_Rho)),
    hjust = 0, size = 2.6, color = "grey25", inherit.aes = FALSE
  ) +
  annotate("text", x = 0.78, y = length(s3_ord6) + 0.65, label = "Fidelity",
           hjust = 0, size = 2.7, fontface = "bold", color = "grey20") +
  scale_color_manual(values = c("0_Unadjusted_Raw" = "#D62728",
                                "3_scRNA_Gene_Purification" = "#00A087"),
                     labels = c("Original (Raw)", "scRNA-Purified (S3)"), name = NULL) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Original (Raw)", "scRNA-Purified (S3)"), name = NULL) +
  coord_cartesian(xlim = c(-0.22, 0.96), clip = "off") +
  labs(title = "Signature-Level Drift Suppression (S3)",
       subtitle = "3-cohort mean at 30% admixture; right = pure-tumor fidelity",
       x = "Mean Standardized Drift (SD)", y = NULL) +
  theme_pub(9.2) +
  theme(legend.position = "bottom")

# Panel E (重磅新增): 789 例真实患者 OS 预后“解封与保真”验证 (Raw vs S3 Purified)
key_unmask <- c("HM_Complement", "HM_IL2_STAT5", "HM_TNFA_via_NFKB", "T_cell_exhaustion",
                "HM_Allograft_Rejection", "TAM_Peng", "HM_IFN_gamma", "HM_IFN_alpha", "MDSC_Peng")

p6e_wide <- res4_cli$meta_cox_df %>%
  select(Signature, HR_Raw, P_Raw) %>%
  inner_join(res5_obj$pur_meta_cox, by = "Signature") %>%
  filter(Signature %in% key_unmask) %>%
  mutate(Sig_Label = map_sig(Signature))

ord_6e <- rev(map_sig(key_unmask))

p6e_long <- bind_rows(
  p6e_wide %>% transmute(Sig_Label, Ver = "Original (Raw)", HR = HR_Raw, P = P_Raw),
  p6e_wide %>% transmute(Sig_Label, Ver = "scRNA-Purified (S3)", HR = HR_Purified_Meta, P = P_Purified_Meta)
) %>%
  mutate(
    Sig_Label = factor(Sig_Label, levels = ord_6e),
    Sig_Stat  = ifelse(P < 0.05, "P < 0.05", "NS"),
    p_lbl     = case_when(P < 0.001 ~ "***", P < 0.01 ~ "**", P < 0.05 ~ "*", TRUE ~ "")
  )

p6e <- ggplot(p6e_long, aes(x = HR, y = Sig_Label)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey65", linewidth = 0.9) +
  geom_point(aes(color = Ver, shape = Sig_Stat), size = 3.0) +
  geom_text(
    data = p6e_long %>% filter(Ver == "scRNA-Purified (S3)", p_lbl != ""),
    aes(label = p_lbl), nudge_y = 0.28, size = 3.2, color = "#00796B", fontface = "bold"
  ) +
  scale_color_manual(values = c("Original (Raw)" = "#D62728",
                                "scRNA-Purified (S3)" = "#00A087"), name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "NS" = 1), name = NULL) +
  labs(title = "Real-World OS Unmasking (n=789)",
       subtitle = "3-cohort Meta-HR; filled = P < 0.05",
       x = "Pooled Hazard Ratio (OS)", y = NULL) +
  theme_pub(9.2) +
  theme(legend.position = "bottom", legend.direction = "vertical")

fig6 <- ((p6a / p6b / p6c) | p6d | p6e) +
  plot_layout(widths = c(0.95, 1.25, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig6, "Figure6_Correction_Benchmark_and_Unmasking", width = 16.0, height = 9.5)

message("\n>>> 全部 6 张投稿级主图精修完成！文件已保存至: ", out_dir)
print(list.files(out_dir, pattern = "\\.png$"))
