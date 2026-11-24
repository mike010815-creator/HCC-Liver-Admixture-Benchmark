# ==============================================================================
# 脚本名称: 06b_finish_fig5_fig6.R
# 核心任务: 修复 patchwork 嵌套布局与 ggplot2 4.0.0 水平误差棒语法，
#           3 秒补齐投稿级 Figure 5 与 Figure 6
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

out_dir <- "04_figures/final_manuscript"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. 与 Figure 1–4 完全一致的规范英文标签、配色与主题
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
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height, device = "pdf")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, dpi = 300, bg = "white")
  tryCatch(
    ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot_obj, width = width, height = height, dpi = 300,
           device = "tiff", compression = "lzw", bg = "white"),
    error = function(e) message("    (TIFF 跳过: ", e$message, ")")
  )
  message(sprintf("    已导出 -> %s (.pdf / .png / .tiff)", stem))
}

res4_cli <- readRDS("01_clean_data/04b_clinical_survival_results.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

# ==============================================================================
# FIGURE 5: 789 例真实初治原发 HCC 共变验证与临床预后偏倚
# ==============================================================================
message(">>> [1/2] 正在生成投稿级 Figure 5...")

real_cor_wide <- res4_cli$real_cor_wide %>%
  mutate(Sig_Label = map_sig(Signature), Axis_Label = map_axis(Axis))
sig_ord5 <- rev(real_cor_wide$Sig_Label)

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
  geom_errorbar(aes(xmin = Lo, xmax = Hi), width = 0.22, linewidth = 0.7,
                orientation = "y", position = position_dodge(width = 0.52)) +
  geom_point(aes(shape = Sig_Flag), size = 2.6, position = position_dodge(width = 0.52)) +
  scale_color_manual(values = c("Low-NHAS (Purer Tumors)"      = "#1F77B4",
                                "High-NHAS (Higher Admixture)" = "#D62728"), name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "NS" = 1), name = NULL) +
  labs(title = "NHAS-Stratified Prognostic Attenuation",
       subtitle = "3-cohort Meta-HR (95% CI); * P_int < 0.05, ** P_int < 0.01",
       x = "Pooled Hazard Ratio (95% CI)", y = NULL) +
  theme_pub(9.5) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# 用 wrap_plots 明确分层，彻底杜绝 "|" 自动展平导致的列数不匹配
left_5  <- wrap_plots(p5a, p5c, ncol = 2, widths = c(1.15, 0.95))
right_5 <- wrap_plots(p5b, p5d, ncol = 1, heights = c(0.95, 1.15))
fig5    <- wrap_plots(left_5, right_5, ncol = 2, widths = c(1.38, 1.12)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig5, "Figure5_Clinical_Prognostic_Bias", width = 15.5, height = 9.5)

# ==============================================================================
# FIGURE 6: 三大纠偏策略基准测试 + 单细胞净化降漂 + 789 例临床预后“解封”验证
# ==============================================================================
message(">>> [2/2] 正在生成投稿级 Figure 6...")

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

p6a <- bench_df %>%
  filter(Scenario == "Uniform_30pct", Signature %in% pruned_sigs) %>%
  ggplot(aes(x = Strategy, y = Abs_Drift, fill = Strategy)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.9, linewidth = 0.45) +
  scale_x_discrete(labels = strat_labels) +
  scale_fill_manual(values = strat_cols, guide = "none") +
  labs(title = "Absolute Drift (Uniform 30%)", x = NULL, y = "|Standardized Drift| (SD)") +
  theme_pub(9.2)

p6b <- bench_df %>%
  filter(Scenario == "Heterogeneous_0_40pct", Signature %in% pruned_sigs) %>%
  ggplot(aes(x = Strategy, y = Spearman_Rho, fill = Strategy)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.9, linewidth = 0.45) +
  scale_x_discrete(labels = strat_labels) +
  scale_fill_manual(values = strat_cols, guide = "none") +
  labs(title = "Rank Stability (Variable 0-40%)", x = NULL, y = "Spearman rho vs. Pure") +
  theme_pub(9.2)

p6c <- bench_df %>%
  filter(Scenario == "Heterogeneous_0_40pct", Signature %in% pruned_sigs) %>%
  ggplot(aes(x = Strategy, y = Flip_Fixed_Pct, fill = Strategy)) +
  geom_boxplot(width = 0.55, alpha = 0.85, outlier.size = 0.9, linewidth = 0.45) +
  scale_x_discrete(labels = strat_labels) +
  scale_fill_manual(values = strat_cols, guide = "none") +
  labs(title = "Misclassification (Variable 0-40%)", x = NULL, y = "Flip Rate (%)") +
  theme_pub(9.2)

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

left_6 <- wrap_plots(p6a, p6b, p6c, ncol = 1)
fig6   <- wrap_plots(left_6, p6d, p6e, ncol = 3, widths = c(0.95, 1.25, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

save_pub_fig(fig6, "Figure6_Correction_Benchmark_and_Unmasking", width = 16.0, height = 9.5)

message("\n>>> 全部 6 张投稿级主图已完整集齐！目录文件清单如下:")
print(list.files(out_dir, pattern = "\\.png$"))
