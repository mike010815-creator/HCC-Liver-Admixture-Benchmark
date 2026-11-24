# ==============================================================================
# 脚本名称: 08_final_figures_4_to_6.R
# 核心任务: 基于 07a/07b 全部补算结果，用统一 R ggplot2 出版体系重绘 Figure 4–6
# 输出目录: 04_figures/final_manuscript_v2/ (.pdf / .png / .tiff)
# ==============================================================================

rm(list = ls())
set.seed(2026)

if (!dir.exists("01_clean_data") && dir.exists("~/HCC_Admixture_Project/01_clean_data")) {
  setwd("~/HCC_Admixture_Project")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(survival)
  library(patchwork)
})

rev_dir <- "03_results/revision_20261123"
out_dir <- "04_figures/final_manuscript_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. 全局出版级主题、规范英文标签与统一配色 (与 Figure 1-3 完全一致)
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
      plot.title        = element_text(face = "bold", size = base_size + 0.8, hjust = 0, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = base_size - 1.3, color = "grey35", margin = margin(b = 5)),
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

pool_meta <- function(beta, se) {
  ok <- !is.na(beta) & !is.na(se) & se > 0 & is.finite(beta) & is.finite(se)
  b <- beta[ok]; s <- se[ok]; k <- length(b)
  w <- 1 / (s^2); B <- sum(w * b) / sum(w); S <- sqrt(1 / sum(w))
  Q <- sum(w * (b - B)^2)
  I2 <- if (k > 1 && Q > 0) max(0, (Q - k + 1) / Q) * 100 else 0
  tau2 <- if (k > 1) max(0, (Q - k + 1) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (s^2 + tau2); Br <- sum(wr * b) / sum(wr); Sr <- sqrt(1 / sum(wr))
  list(
    HR = exp(B), Lo = exp(B - 1.96 * S), Hi = exp(B + 1.96 * S), P = 2 * pnorm(-abs(B / S)),
    I2 = I2,
    R_HR = exp(Br), R_Lo = exp(Br - 1.96 * Sr), R_Hi = exp(Br + 1.96 * Sr), R_P = 2 * pnorm(-abs(Br / Sr))
  )
}

# 加载基础对象
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
res4_cli <- readRDS("01_clean_data/04b_clinical_survival_results.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

sig_list     <- sigs_obj$sig_list
sig_meta     <- sigs_obj$sig_meta
nhas_genes   <- res3_obj$nhas_genes
pur_sig_list <- res5_obj$purified_sig_list
pruned_sigs  <- res5_obj$prune_summary$Signature[res5_obj$prune_summary$Was_Pruned]

# ==============================================================================
# FIGURE 4: 全 7 张 Visium 空间切片统一分析 + 代表性切片 HCC3R 空间分布
# ==============================================================================
message(">>> [1/3] 正在生成升级版 Figure 4 (全部 7 张空间切片统一流程)...")

sp7_qc  <- read.csv(file.path(rev_dir, "07b_spatial_all7_slides_qc.csv"))
sp7_cor <- read.csv(file.path(rev_dir, "07b_spatial_all7_slides_cor.csv"))
sp7_sum <- read.csv(file.path(rev_dir, "07b_spatial_all7_sig_summary.csv")) %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature") %>%
  mutate(Sig_Label = map_sig(Signature), Axis_Label = map_axis(Axis))

# 读取 HCC3R 生成 4 张代表性空间彩图 (A–D)
h3_dir <- "01_clean_data/spatial_unpacked/GSE238264/HCC3R"
stopifnot(dir.exists(h3_dir))
cnt3 <- Read10X(file.path(h3_dir, "filtered_feature_bc_matrix"))
if (is.list(cnt3)) cnt3 <- cnt3[[1]]
seu3 <- CreateSeuratObject(counts = cnt3, assay = "Spatial", project = "HCC3R")
img3 <- Read10X_Image(file.path(h3_dir, "spatial"), image.name = "tissue_lowres_image.png", filter.matrix = TRUE)
DefaultAssay(img3) <- "Spatial"
bc3 <- intersect(colnames(seu3), Cells(img3))
seu3 <- subset(seu3, cells = bc3)
img3 <- img3[bc3]
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
seu3$SP_Comp_Pur <- score_ms(seu3, pur_sig_list$HM_Complement, "M_COMPP_")

sp_specs <- c(
  "SP_NHAS"     = "NHAS (HCC3R)",
  "SP_G2M"      = "G2M Checkpoint (Prolif. Ctrl)",
  "SP_Comp"     = "Hallmark Complement (Raw)",
  "SP_Comp_Pur" = "Hallmark Complement (S3 Pruned)"
)

sp_maps <- lapply(names(sp_specs), function(cn) {
  med_v <- median(seu3[[cn]][, 1], na.rm = TRUE)
  suppressMessages(
    SpatialFeaturePlot(seu3, features = cn, image.alpha = 0.35, pt.size.factor = 1.6) +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                           midpoint = med_v, name = "Score") +
      ggtitle(sp_specs[[cn]]) +
      theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0.5),
            legend.position = "bottom",
            legend.key.width = unit(0.65, "cm"),
            legend.key.height = unit(0.28, "cm"),
            legend.title = element_text(size = 7.8, face = "bold"),
            legend.text = element_text(size = 7.2),
            plot.margin = margin(2, 4, 2, 4))
  )
})
p4_top <- wrap_plots(sp_maps, ncol = 4)

# Panel E: 全部 7 张切片 × 24 特征统一热图 (切片按 NHAS span 从宽到窄排序)
slide_ord <- sp7_qc$Slide
sig_ord4  <- rev(sp7_sum$Sig_Label)

p4e <- sp7_cor %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = sig_ord4),
    Slide     = factor(Slide, levels = slide_ord)
  ) %>%
  ggplot(aes(x = Slide, y = Sig_Label, fill = Spatial_Rho_Raw)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Spatial_Rho_Raw)), size = 2.55, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-0.45, 0.60), oob = scales::squish,
                       name = "Spatial rho\nwith NHAS") +
  labs(
    title    = "All 7 Visium Sections: Background-Adjusted NHAS Co-variation",
    subtitle = "Unified AddModuleScore pipeline (spots >= 500 features); slides ordered by NHAS span",
    x = NULL, y = NULL
  ) +
  theme_pub(9.2) +
  theme(axis.text.x = element_text(face = "bold", angle = 20, hjust = 1),
        legend.position = "bottom", legend.key.width = unit(1.0, "cm"))

# Panel F: 7 切片中位空间相关 (含 IQR) vs 三队列 Bulk 30% 漂移一致性散点图
conc_7 <- cor(sp7_sum$Median_Rho_7slides, sp7_sum$Bulk30_Mean, method = "spearman")
sp7_iqr <- sp7_cor %>%
  group_by(Signature) %>%
  summarise(Q25 = quantile(Spatial_Rho_Raw, 0.25, na.rm = TRUE),
            Q75 = quantile(Spatial_Rho_Raw, 0.75, na.rm = TRUE), .groups = "drop")
p4f_df <- sp7_sum %>% left_join(sp7_iqr, by = "Signature")

label_sigs4 <- c("HM_Complement", "CTRL_G2M_Checkpoint", "T_cell_exhaustion",
                 "HM_IL6_JAK_STAT3", "TAM_Peng", "HM_EMT", "Treg_Signature")
p4f_lbl <- p4f_df %>%
  filter(Signature %in% label_sigs4) %>%
  mutate(Short = recode(Signature,
                        "HM_Complement"       = "Complement (7/7 pos.)",
                        "CTRL_G2M_Checkpoint" = "G2M Ctrl (0/7 pos.)",
                        "T_cell_exhaustion"   = "T-cell Exhaustion",
                        "HM_IL6_JAK_STAT3"    = "IL6-JAK-STAT3",
                        "TAM_Peng"            = "TAM",
                        "HM_EMT"              = "EMT",
                        "Treg_Signature"      = "Treg"))

p4f <- ggplot(p4f_df, aes(x = Bulk30_Mean, y = Median_Rho_7slides, color = Axis_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.75) +
  geom_errorbar(aes(ymin = Q25, ymax = Q75), width = 0.015, alpha = 0.45, linewidth = 0.45) +
  geom_point(size = 2.8, alpha = 0.92) +
  geom_text(
    data = p4f_lbl,
    aes(label = Short),
    size = 2.75, fontface = "bold", color = "grey15",
    nudge_y = ifelse(p4f_lbl$Signature %in% c("HM_Complement", "TAM_Peng", "HM_IL6_JAK_STAT3"), 0.025, -0.028),
    hjust = ifelse(p4f_lbl$Bulk30_Mean > 0.55, 1.0, 0.1),
    check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_color_manual(values = axis_cols, name = "Biological Axis") +
  annotate("text", x = -0.36, y = 0.36, hjust = 0, vjust = 1, fontface = "bold", size = 3.2,
           label = sprintf("7-Slide Median vs. Bulk 30%% Drift\nSpearman rho = %.3f (24 gene sets)\nComplement: 7/7 positive | G2M: 7/7 negative",
                           conc_7)) +
  coord_cartesian(ylim = c(-0.28, 0.39), xlim = c(-0.40, 0.76)) +
  labs(
    title    = "Cross-Modal Concordance Across All 7 Sections",
    subtitle = "Points = median across 7 slides; vertical bars = interquartile range (IQR)",
    x        = "Mean 3-Cohort Bulk 30% Admixture Drift (SD)",
    y        = "7-Slide Median Spatial rho with NHAS"
  ) +
  theme_pub(9.2) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom")

p4_bot <- wrap_plots(p4e, p4f, ncol = 2, widths = c(1.18, 1.08))
fig4 <- wrap_plots(p4_top, p4_bot, ncol = 1, heights = c(0.78, 1.55)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12.5))

save_pub_fig(fig4, "Figure4_Spatial_All7_Validation", width = 15.0, height = 11.8)

# ==============================================================================
# FIGURE 5: 789 例真实初治原发 HCC 共变、统一 SD 分层与基质四分位非线性诊断
# ==============================================================================
message(">>> [2/3] 正在生成升级版 Figure 5 (全队列统一 SD 分层 + 基质 Q1-Q4 非线性诊断)...")

real_cor_wide <- res4_cli$real_cor_wide %>%
  mutate(Sig_Label = map_sig(Signature), Axis_Label = map_axis(Axis))
sig_ord5 <- rev(real_cor_wide$Sig_Label)

# Panel A: 三队列真实肿瘤 NHAS 相关热图
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
  geom_text(aes(label = sprintf("%.2f%s", Real_Rho_with_NHAS, star)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Spearman rho\nwith NHAS") +
  labs(title = "Real-Tumor Co-variation with NHAS",
       subtitle = "789 primary tumors (nominal P stars)", x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(axis.text.x = element_text(face = "bold", size = 8.2),
        legend.position = "bottom", legend.key.width = unit(0.95, "cm"))

# Panel B: Raw vs NHAS-Adjusted Meta-Cox OS HR
meta_cox <- res4_cli$meta_cox_df %>% mutate(Sig_Label = map_sig(Signature))
p5b <- meta_cox %>%
  select(Sig_Label, HR_Raw, HR_Adj) %>%
  pivot_longer(c(HR_Raw, HR_Adj), names_to = "Model", values_to = "HR") %>%
  mutate(
    Sig_Label = factor(Sig_Label, levels = sig_ord5),
    Model = recode(Model, HR_Raw = "Unadjusted (~ Sig)",
                   HR_Adj = "NHAS-Adjusted (~ Sig + NHAS)")
  ) %>%
  ggplot(aes(x = HR, y = Sig_Label)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_line(aes(group = Sig_Label), color = "grey65", linewidth = 0.85) +
  geom_point(aes(color = Model, shape = Model), size = 2.65) +
  scale_color_manual(values = c("Unadjusted (~ Sig)" = "#D62728",
                                "NHAS-Adjusted (~ Sig + NHAS)" = "#1F77B4"), name = NULL) +
  scale_shape_manual(values = c(17, 16), name = NULL) +
  labs(title = "NHAS Adjustment Shift (n=789)",
       subtitle = "Fixed-effect Meta-HR per 1-SD",
       x = "Pooled Hazard Ratio (OS)", y = NULL) +
  theme_pub(9.2) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "vertical")

# Panel C: 全队列统一 SD 分层 Meta-Cox (现算现用，确保量纲严格一致)
key_immune_strat <- c("CYT_Rooney", "CD8_T_effector", "T_cell_inflamed_GEP",
                      "HM_IFN_gamma", "HM_IFN_alpha", "Immune_Checkpoint",
                      "TAM_Peng", "MDSC_Peng", "CTRL_G2M_Checkpoint")

common_sd_cohort <- map_dfr(names(res4_cli$real_all_list), function(cn) {
  d <- res4_cli$real_all_list[[cn]]
  med_n <- median(d$NHAS)
  map_dfr(key_immune_strat, function(s) {
    d$Sig_z <- as.numeric(scale(d[[s]]))  # 全队列统一 1-SD
    map_dfr(c("Low NHAS", "High NHAS"), function(st) {
      sub <- if (st == "Low NHAS") d[d$NHAS < med_n, ] else d[d$NHAS >= med_n, ]
      sm  <- summary(coxph(Surv(OS_time, OS_event) ~ Sig_z, data = sub))$coefficients
      tibble(Cohort = cn, Signature = s, Stratum = st,
             beta = unname(sm[1, "coef"]), se = unname(sm[1, "se(coef)"]))
    })
  })
})

common_sd_meta <- common_sd_cohort %>%
  group_by(Signature, Stratum) %>%
  summarise(res = list(pool_meta(beta, se)), .groups = "drop") %>%
  unnest_wider(res) %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = rev(map_sig(key_immune_strat))),
    Stratum   = factor(Stratum, levels = c("Low NHAS", "High NHAS")),
    Sig_Flag  = ifelse(P < 0.05, "P < 0.05", "NS")
  )

p5c <- ggplot(common_sd_meta, aes(x = HR, y = Sig_Label, color = Stratum)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_errorbar(aes(xmin = Lo, xmax = Hi), width = 0.22, linewidth = 0.65,
                orientation = "y", position = position_dodge(width = 0.52)) +
  geom_point(aes(shape = Sig_Flag), size = 2.5, position = position_dodge(width = 0.52)) +
  scale_color_manual(values = c("Low NHAS" = "#1F77B4", "High NHAS" = "#D62728"), name = NULL) +
  scale_shape_manual(values = c("P < 0.05" = 16, "NS" = 1), name = NULL) +
  labs(
    title    = "NHAS-Stratified Associations (Common Cohort SD)",
    subtitle = "Exploratory 3-cohort Meta-HR (95% CI); overlapping CIs indicate no confirmed subgroup shift",
    x        = "Pooled Hazard Ratio (95% CI)", y = NULL
  ) +
  theme_pub(9.2) +
  theme(legend.position = "bottom")

# Panel D: 基质通路 NHAS 四分位 (Q1->Q4) 非单调轨迹图 (揭示左尾极低分化杠杆效应)
q_meta <- read.csv(file.path(rev_dir, "07b_stromal_quartile_HRs.csv")) %>%
  mutate(
    Sig_Label = factor(map_sig(Signature), levels = c("Pan-F-TBRs (TGF-beta)", "Hallmark EMT", "CAF")),
    NHAS_Q    = factor(NHAS_Q, levels = c("Q1_Lowest", "Q2", "Q3", "Q4_Highest"),
                       labels = c("Q1\n(Lowest)", "Q2", "Q3", "Q4\n(Highest)"))
  )

p5d <- ggplot(q_meta, aes(x = NHAS_Q, y = HR, group = Sig_Label, color = Sig_Label)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_errorbar(aes(ymin = Lo, ymax = Hi), width = 0.14, linewidth = 0.6,
                position = position_dodge(width = 0.32)) +
  geom_line(linewidth = 0.85, position = position_dodge(width = 0.32)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.32)) +
  scale_color_manual(values = c("Pan-F-TBRs (TGF-beta)" = "#3C5488",
                                "Hallmark EMT"          = "#E64B35",
                                "CAF"                   = "#00A087"), name = NULL) +
  labs(
    title    = "Non-Monotonic Stromal HR Across NHAS Quartiles",
    subtitle = "Left-skewed NHAS (skew -2.3 to -2.5) drives apparent linear interaction; binary interaction NS (P > 0.24)",
    x        = "NHAS Expression Quartile (Full-Cohort SD)",
    y        = "Pooled OS Hazard Ratio (95% CI)"
  ) +
  theme_pub(9.2) +
  theme(legend.position = "bottom")

left_5  <- wrap_plots(p5a, p5b, ncol = 2, widths = c(1.18, 0.92))
right_5 <- wrap_plots(p5c, p5d, ncol = 1, heights = c(1.15, 0.95))
fig5    <- wrap_plots(left_5, right_5, ncol = 2, widths = c(1.30, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12.5))

save_pub_fig(fig5, "Figure5_Clinical_Survival_and_NHAS_Diagnostics", width = 15.8, height = 9.8)

# ==============================================================================
# FIGURE 6: 基准对比 + 表达匹配随机对照 + Bootstrap 分层权衡 + 多变量生存森林图
# ==============================================================================
message(">>> [3/3] 正在生成升级版 Figure 6 (含表达匹配随机对照、Bootstrap 95% CI 与多变量森林图)...")

boot_df  <- read.csv(file.path(rev_dir, "07a_patient_bootstrap_95CI.csv"))
null_g   <- read.csv(file.path(rev_dir, "07a_random_deletion_null_grand.csv"))
null_sig <- read.csv(file.path(rev_dir, "07a_random_deletion_null_bysig.csv"))
ext_s2   <- read.csv(file.path(rev_dir, "07a_external_S2_validation.csv"))
mv_coh   <- read.csv(file.path(rev_dir, "07b_multivar_pruned_cohort.csv"))

bench_u30 <- res5_obj$bench_df %>%
  filter(Scenario == "Uniform_30pct", Signature %in% pruned_sigs) %>%
  group_by(Strategy) %>%
  summarise(Mean_Drift = mean(Abs_Drift), Mean_Rho = mean(Spearman_Rho),
            Mean_FFlip = mean(Flip_Fixed_Pct), .groups = "drop")

b_all <- boot_df %>% filter(Group == "All_18_Pruned")
ng_ex <- null_g %>% filter(Null_Type == "ExprBin_SizeMatched")

# Panel A: 均匀 30% 混入横向对比 (Raw / S1 / S2 内 / S2 外 / 表达匹配随机对照 / S3)
p6a_df <- tibble(
  Method = c("Raw\n(Unadj.)", "S1: NHAS\nResidual", "S2: Linear\n(Within)",
             "S2: Linear\n(Ext. I->T)", "Random Null\n(Expr-Match)", "S3: scRNA\nPruning"),
  Val    = c(
    b_all$Estimate[b_all$Metric == "raw_abs_drift"],
    bench_u30$Mean_Drift[bench_u30$Strategy == "1_NHAS_Residualization"],
    bench_u30$Mean_Drift[bench_u30$Strategy == "2_Linear_Ref_Subtraction"],
    ext_s2$Abs_Drift[ext_s2$Direction == "ICGC_LIRI_JP -> TCGA_LIHC" & ext_s2$Method == "External_S2"],
    ng_ex$Null_AbsDrift_Mean,
    b_all$Estimate[b_all$Metric == "s3_abs_drift"]
  ),
  Lo = c(b_all$CI_Low[b_all$Metric == "raw_abs_drift"], NA, NA, NA,
         ng_ex$Null_AbsDrift_P05, b_all$CI_Low[b_all$Metric == "s3_abs_drift"]),
  Hi = c(b_all$CI_High[b_all$Metric == "raw_abs_drift"], NA, NA, NA,
         ng_ex$Null_AbsDrift_P95, b_all$CI_High[b_all$Metric == "s3_abs_drift"]),
  Fill = c("#D62728", "#FF7F0E", "#1F77B4", "#6BAED6", "#7F7F7F", "#00A087")
) %>%
  mutate(Method = factor(Method, levels = Method))

p6a <- ggplot(p6a_df, aes(x = Method, y = Val, fill = Method)) +
  geom_col(width = 0.62, alpha = 0.88) +
  geom_errorbar(aes(ymin = Lo, ymax = Hi), width = 0.2, linewidth = 0.6, na.rm = TRUE) +
  geom_text(aes(y = ifelse(is.na(Hi), Val, Hi) + 0.05, label = sprintf("%.3f", Val)),
            size = 2.8, fontface = "bold") +
  annotate("text", x = 5.5, y = 1.24, hjust = 0.5, size = 2.85, fontface = "bold", color = "#00796B",
           label = "S3 vs. Matched Null:\nP_emp = 0.0196 (50/50)\nFlip: 15.0% vs. 20.6%") +
  scale_fill_manual(values = setNames(p6a_df$Fill, p6a_df$Method), guide = "none") +
  coord_cartesian(ylim = c(0, 1.50)) +
  labs(
    title    = "Uniform 30% Admixture: Strategy & Null Control Benchmark",
    subtitle = "Bars = 18-signature mean |Drift|; error bars = 95% bootstrap CI (Raw/S3) or 5-95% null range",
    x = NULL, y = "Mean |Standardized Drift| (SD)"
  ) +
  theme_pub(9.0) +
  theme(axis.text.x = element_text(size = 8.0))

# Panel B: 患者配对 Bootstrap 95% CI —— 自身稳健性 vs 共同原始目标保真度 (按基线一致度分层)
p6b_df <- boot_df %>%
  filter(Group %in% c("All_18_Pruned", "High_Fidelity_ge09", "Redesigned_lt09"),
         Metric %in% c("raw_rho", "s3_self_rho", "s3_orig_rho")) %>%
  mutate(
    Group_lbl = factor(Group,
                       levels = c("All_18_Pruned", "High_Fidelity_ge09", "Redesigned_lt09"),
                       labels = c("All 18\nPruned", "High-Fidelity\n(n=11, rho>=0.90)", "Redesigned\n(n=7, rho<0.90)")),
    Target_lbl = factor(Metric,
                        levels = c("raw_rho", "s3_self_rho", "s3_orig_rho"),
                        labels = c("Raw (vs. Original Baseline)",
                                   "S3 (vs. Own Pruned Baseline)",
                                   "S3 (vs. Original Baseline)"))
  )

p6b <- ggplot(p6b_df, aes(x = Group_lbl, y = Estimate, color = Target_lbl, shape = Target_lbl)) +
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.16, linewidth = 0.7,
                position = position_dodge(width = 0.48)) +
  geom_point(size = 2.9, position = position_dodge(width = 0.48)) +
  geom_text(aes(label = sprintf("%.3f", Estimate), y = CI_High + 0.015),
            size = 2.55, position = position_dodge(width = 0.48), show.legend = FALSE) +
  scale_color_manual(values = c("Raw (vs. Original Baseline)"  = "#D62728",
                                "S3 (vs. Own Pruned Baseline)" = "#00A087",
                                "S3 (vs. Original Baseline)"   = "#3C5488"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  coord_cartesian(ylim = c(0.68, 0.985)) +
  labs(
    title    = "Self-Stability vs. Original-Target Fidelity Trade-off",
    subtitle = "Patient-paired bootstrap 95% CI (B = 1,000, uniform 30% admixture)",
    x = NULL, y = "Rank Correlation (Spearman rho)"
  ) +
  theme_pub(9.0) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# Panel C: 18 个特征逐一降漂图 (Raw vs 表达匹配随机对照 vs S3, 附单特征 P_emp 星号与基线 Fidelity)
raw_by_sig <- res5_obj$bench_df %>%
  filter(Scenario == "Uniform_30pct", Strategy == "0_Unadjusted_Raw", Signature %in% pruned_sigs) %>%
  group_by(Signature) %>%
  summarise(Raw_Abs = mean(Abs_Drift), .groups = "drop")

p6c_wide <- null_sig %>%
  left_join(raw_by_sig, by = "Signature") %>%
  left_join(res5_obj$fidelity_df, by = "Signature") %>%
  mutate(
    Sig_Label = map_sig(Signature),
    Sig_Star  = ifelse(Emp_P_Drift < 0.05, paste0(Sig_Label, " *"), Sig_Label),
    Gain_vs_Null = Null_Drift - S3_Drift
  ) %>%
  arrange(desc(Gain_vs_Null))

ord_6c <- rev(p6c_wide$Sig_Star)

p6c_long <- p6c_wide %>%
  select(Sig_Star, Raw_Abs, Null_Drift, S3_Drift) %>%
  pivot_longer(-Sig_Star, names_to = "Type", values_to = "Abs_Drift") %>%
  mutate(
    Sig_Star = factor(Sig_Star, levels = ord_6c),
    Type = factor(Type,
                  levels = c("Raw_Abs", "Null_Drift", "S3_Drift"),
                  labels = c("Original (Raw)", "Expr-Matched Random Null (R=50)", "S3 scRNA-Pruned"))
  )

p6c <- ggplot(p6c_long, aes(x = Abs_Drift, y = Sig_Star)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_segment(
    data = p6c_wide %>% mutate(Sig_Star = factor(Sig_Star, levels = ord_6c)),
    aes(x = S3_Drift, xend = Raw_Abs, y = Sig_Star, yend = Sig_Star),
    color = "grey70", linewidth = 0.85, inherit.aes = FALSE
  ) +
  geom_point(aes(color = Type, shape = Type), size = 2.75) +
  geom_text(
    data = p6c_wide %>% mutate(Sig_Star = factor(Sig_Star, levels = ord_6c)),
    aes(x = 0.77, y = Sig_Star, label = sprintf("%.2f", Mean_PureTumor_Fidelity_Rho)),
    hjust = 0, size = 2.6, color = "grey25", inherit.aes = FALSE
  ) +
  annotate("text", x = 0.77, y = length(ord_6c) + 0.90, label = "Base rho",
           hjust = 0, vjust = 0.5, size = 2.75, fontface = "bold", color = "grey15") +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.45))) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6)) +
  scale_color_manual(values = c("Original (Raw)"                 = "#D62728",
                                "Expr-Matched Random Null (R=50)" = "#7F7F7F",
                                "S3 scRNA-Pruned"                 = "#00A087"), name = NULL) +
  scale_shape_manual(values = c(16, 18, 17), name = NULL) +
  coord_cartesian(xlim = c(0, 0.92)) +
  labs(
    title    = "Signature-Level Drift vs. Expression-Matched Null",
    subtitle = "* P_emp < 0.05 vs. matched random deletion (11/18 signatures); right = baseline rho",
    x        = "3-Cohort Mean |Standardized Drift| at 30% (SD)", y = NULL
  ) +
  theme_pub(9.0) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# Panel D: 763 例完整病例多变量 Cox 敏感性森林图 (调整 Age + Sex + TNM Stage)
mv_all18 <- mv_coh %>%
  group_by(Signature) %>%
  summarise(
    raw_m = list(pool_meta(beta_raw_mv, se_raw_mv)),
    pur_m = list(pool_meta(beta_pur_mv, se_pur_mv)),
    .groups = "drop"
  ) %>%
  mutate(
    Raw_HR = sapply(raw_m, `[[`, "HR"), Raw_Lo = sapply(raw_m, `[[`, "Lo"), Raw_Hi = sapply(raw_m, `[[`, "Hi"),
    Pur_HR = sapply(pur_m, `[[`, "HR"), Pur_Lo = sapply(pur_m, `[[`, "Lo"), Pur_Hi = sapply(pur_m, `[[`, "Hi"),
    Pur_P  = sapply(pur_m, `[[`, "P"),  Pur_I2 = sapply(pur_m, `[[`, "I2"),
    Pur_RHR = sapply(pur_m, `[[`, "R_HR"), Pur_RLo = sapply(pur_m, `[[`, "R_Lo"), Pur_RHi = sapply(pur_m, `[[`, "R_Hi"),
    Pur_RP  = sapply(pur_m, `[[`, "R_P")
  ) %>%
  select(-raw_m, -pur_m) %>%
  mutate(Pur_BH_q = p.adjust(Pur_P, "BH"))

key_mv <- c("HM_Complement", "HM_IL2_STAT5", "T_cell_exhaustion",
            "MDSC_Peng", "HM_IFN_gamma", "HM_IFN_alpha", "TAM_Peng")

p6d_sub <- mv_all18 %>%
  filter(Signature %in% key_mv) %>%
  mutate(
    Sig_Label = map_sig(Signature),
    Anno_Txt  = sprintf("q=%.3f, I2=%.0f%%", Pur_BH_q, Pur_I2)
  )

ord_6d <- rev(map_sig(key_mv))

p6d_long <- bind_rows(
  p6d_sub %>% transmute(Sig_Label, Model = "Raw MV (Fixed-Effect)",       HR = Raw_HR,  Lo = Raw_Lo,  Hi = Raw_Hi),
  p6d_sub %>% transmute(Sig_Label, Model = "S3 Pruned MV (Fixed-Effect)", HR = Pur_HR,  Lo = Pur_Lo,  Hi = Pur_Hi),
  p6d_sub %>% transmute(Sig_Label, Model = "S3 Pruned MV (Random-Eff.)",  HR = Pur_RHR, Lo = Pur_RLo, Hi = Pur_RHi)
) %>%
  mutate(
    Sig_Label = factor(Sig_Label, levels = ord_6d),
    Model     = factor(Model, levels = c("Raw MV (Fixed-Effect)",
                                         "S3 Pruned MV (Fixed-Effect)",
                                         "S3 Pruned MV (Random-Eff.)"))
  )

p6d <- ggplot(p6d_long, aes(x = HR, y = Sig_Label, color = Model, shape = Model)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_errorbar(aes(xmin = Lo, xmax = Hi), width = 0.18, linewidth = 0.65,
                orientation = "y", position = position_dodge(width = 0.56)) +
  geom_point(size = 2.55, position = position_dodge(width = 0.56)) +
  geom_text(
    data = p6d_sub %>% mutate(Sig_Label = factor(Sig_Label, levels = ord_6d)),
    aes(x = 1.62, y = Sig_Label, label = Anno_Txt),
    hjust = 0, size = 2.55, color = "grey25", inherit.aes = FALSE
  ) +
  annotate("text", x = 1.62, y = length(ord_6d) + 0.62, label = "S3 Fixed q / I2",
           hjust = 0, vjust = 0.5, size = 2.65, fontface = "bold", color = "grey15") +
  scale_y_discrete(expand = expansion(add = c(0.55, 1.05))) +
  scale_x_continuous(breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)) +
  scale_color_manual(values = c("Raw MV (Fixed-Effect)"       = "#D62728",
                                "S3 Pruned MV (Fixed-Effect)" = "#00A087",
                                "S3 Pruned MV (Random-Eff.)"  = "#3C5488"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  coord_cartesian(xlim = c(0.52, 2.02)) +
  labs(
    title    = "Multivariable Clinical OS Sensitivity (n=763)",
    subtitle = "Adjusted for Age, Sex, and TNM Stage (III-IV vs I-II); 95% CI",
    x        = "Pooled Adjusted Hazard Ratio (OS per 1-SD)", y = NULL
  ) +
  theme_pub(9.0) +
  theme(legend.position = "bottom", legend.direction = "vertical")

top_6 <- wrap_plots(p6a, p6b, ncol = 2, widths = c(1.22, 1.0))
bot_6 <- wrap_plots(p6c, p6d, ncol = 2, widths = c(1.12, 1.08))
fig6  <- wrap_plots(top_6, bot_6, ncol = 1, heights = c(0.88, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12.5))

save_pub_fig(fig6, "Figure6_Benchmark_Null_and_Multivar_Survival", width = 15.5, height = 12.2)

message("\n>>> 全部 3 张升级主图已生成至 ", out_dir)
print(list.files(out_dir, pattern = "\\.png$"))
