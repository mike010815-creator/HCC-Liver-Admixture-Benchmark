# ==============================================================================
# 脚本名称: 06c_fix_fig6.R
# 核心任务: 精修 Figure 6 Panel D 的 Fidelity 标题位置，消除与顶框线重叠，
#           直接覆盖更新 04_figures/final_manuscript/Figure6_*
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

out_dir <- "04_figures/final_manuscript"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

map_sig <- function(x) ifelse(x %in% names(pretty_sig), unname(pretty_sig[x]), x)

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

res4_cli <- readRDS("01_clean_data/04b_clinical_survival_results.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

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

# --- Panel D (重点精修版) ---
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
    aes(x = 0.77, y = Sig_Label, label = sprintf("rho=%.2f", Mean_PureTumor_Fidelity_Rho)),
    hjust = 0, size = 2.65, color = "grey25", inherit.aes = FALSE
  ) +
  # 将 Fidelity 稳妥置于框内顶部留白区，远离顶边框线
  annotate("text", x = 0.77, y = length(s3_ord6) + 0.90, label = "Fidelity",
           hjust = 0, vjust = 0.5, size = 2.85, fontface = "bold", color = "grey15") +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.45))) +
  scale_x_continuous(breaks = c(-0.2, 0.0, 0.2, 0.4, 0.6)) +
  scale_color_manual(values = c("0_Unadjusted_Raw" = "#D62728",
                                "3_scRNA_Gene_Purification" = "#00A087"),
                     labels = c("Original (Raw)", "scRNA-Purified (S3)"), name = NULL) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Original (Raw)", "scRNA-Purified (S3)"), name = NULL) +
  coord_cartesian(xlim = c(-0.22, 0.98)) +
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

stem <- "Figure6_Correction_Benchmark_and_Unmasking"
ggsave(file.path(out_dir, paste0(stem, ".pdf")),  fig6, width = 16.0, height = 9.5, device = "pdf")
ggsave(file.path(out_dir, paste0(stem, ".png")),  fig6, width = 16.0, height = 9.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(stem, ".tiff")), fig6, width = 16.0, height = 9.5, dpi = 300,
       device = "tiff", compression = "lzw", bg = "white")

message(">>> Figure 6 已重新精修并覆盖保存完毕！")
