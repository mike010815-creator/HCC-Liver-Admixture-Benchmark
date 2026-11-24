# ==============================================================================
# 脚本名称: 03b_scrna_physiological_calibration.R
# 核心任务:
#   1. 校正 scRNA-seq 酶解肝细胞丢失偏倚，构建三臂对照:
#      - Arm 1: 纯正常肝细胞 (100% Hep)
#      - Arm 2: 生理比例重构癌旁肝 (80% Hep + 20% Non-Parenchymal)
#      - Arm 3: 纯癌旁非实质微环境 (0% Hep, 100% 免疫/内皮/成纤维)
#   2. 连续扫描肝细胞占比 (0% -> 100%)，检验与真实 TCGA/ICGC Bulk 漂移图谱的匹配度
#   3. 核查 CYT_Rooney / T_cell_exhaustion 在纯肝细胞中呈正漂移的具体原因
#   4. 预览 GSE238264 空间压缩包内部目录，为第 4 关铺路
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(GSVA)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 秒速读取已缓存的稀疏单细胞数据、免疫特征库与 Bulk 结果
# ------------------------------------------------------------------------------
sc_obj    <- readRDS("01_clean_data/gse149614_sparse.rds")
sigs_obj  <- readRDS("01_clean_data/immune_signatures.rds")
bulk_sim  <- readRDS("01_clean_data/02d_multicohort_sim_results.rds")
res3_prev <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")

counts_sp  <- sc_obj$counts
meta       <- sc_obj$meta
sig_list   <- sigs_obj$sig_list
paired_pts <- res3_prev$paired_pts

is_T   <- meta$site == "Tumor"
is_N   <- meta$site == "Normal"
is_Hep <- meta$celltype == "Hepatocyte"

to_cpm <- function(cell_ids) {
  v <- Matrix::rowSums(counts_sp[, cell_ids, drop = FALSE])
  v / sum(v) * 1e6
}

P_T      <- sapply(paired_pts, function(p) to_cpm(meta$cell_id[meta$patient == p & is_T]))
P_Hep    <- sapply(paired_pts, function(p) to_cpm(meta$cell_id[meta$patient == p & is_N & is_Hep]))
P_NonHep <- sapply(paired_pts, function(p) to_cpm(meta$cell_id[meta$patient == p & is_N & !is_Hep]))

# 计算 ssGSEA (纯单样本) 与 Zscore 辅助函数
run_ss <- function(mat_log2) {
  p <- ssgseaParam(as.matrix(mat_log2), sig_list, minSize = 2, normalize = FALSE)
  gsva(p, verbose = FALSE)
}

base_log <- log2(P_T + 1)
s0_ss    <- run_ss(base_log)

eval_n_ref <- function(N_ref, label, alpha = 0.30) {
  m_mix <- log2((1 - alpha) * P_T + alpha * N_ref + 1)
  sa_ss <- run_ss(m_mix)
  map_dfr(rownames(s0_ss), function(sig) {
    v0 <- as.numeric(s0_ss[sig, ]); va <- as.numeric(sa_ss[sig, ])
    tibble(
      Arm       = label,
      Signature = sig,
      Std_Drift = mean((va - v0) / sd(v0))
    )
  })
}

# ------------------------------------------------------------------------------
# 2. 扫描肝细胞比例 f_hep (0% -> 100%)，检验与真实 TCGA / ICGC Bulk 漂移的匹配度
# ------------------------------------------------------------------------------
bulk_tcga <- bulk_sim %>% filter(Method == "ssGSEA", Alpha_Pct == 30, grepl("TCGA", Cohort))
bulk_icgc <- bulk_sim %>% filter(Method == "ssGSEA", Alpha_Pct == 30, grepl("ICGC", Cohort))

hep_fracs <- seq(0, 1, by = 0.1)
message(">>> [1/4] 正在扫描癌旁肝细胞比例 (0% -> 100%) 与真实 Bulk 队列的一致性...")

sweep_list <- lapply(hep_fracs, function(fh) {
  N_cal <- fh * P_Hep + (1 - fh) * P_NonHep
  eval_n_ref(N_cal, paste0("Hep_", fh * 100)) %>% mutate(Hep_Frac = fh * 100)
})
sweep_df <- bind_rows(sweep_list)

sweep_cor <- sweep_df %>%
  left_join(bulk_tcga %>% select(Signature, TCGA_Drift = Std_Drift), by = "Signature") %>%
  left_join(bulk_icgc %>% select(Signature, ICGC_Drift = Std_Drift), by = "Signature") %>%
  group_by(Hep_Frac) %>%
  summarise(
    Rho_vs_TCGA = cor(Std_Drift, TCGA_Drift, method = "spearman"),
    Rho_vs_ICGC = cor(Std_Drift, ICGC_Drift, method = "spearman"),
    R_vs_TCGA   = cor(Std_Drift, TCGA_Drift, method = "pearson"),
    R_vs_ICGC   = cor(Std_Drift, ICGC_Drift, method = "pearson"),
    .groups     = "drop"
  )

print(as.data.frame(sweep_cor %>% mutate(across(-Hep_Frac, ~ round(.x, 3)))), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. 提取三臂对比结果: 100% 肝细胞 vs 80% 生理重构癌旁 vs 0% 纯非实质微环境
# ------------------------------------------------------------------------------
tri_df <- sweep_df %>%
  filter(Hep_Frac %in% c(100, 80, 0)) %>%
  mutate(Arm_Label = case_when(
    Hep_Frac == 100 ~ "1_Pure_Hepatocytes (100% Hep)",
    Hep_Frac == 80  ~ "2_Physiological_Liver (80% Hep + 20% NonHep)",
    Hep_Frac == 0   ~ "3_Pure_NonParenchymal (0% Hep, Immune/Stroma)"
  ))

tri_wide <- tri_df %>%
  select(Signature, Hep_Frac, Std_Drift) %>%
  pivot_wider(names_from = Hep_Frac, values_from = Std_Drift, names_prefix = "Hep_") %>%
  left_join(bulk_tcga %>% select(Signature, TCGA_Bulk = Std_Drift), by = "Signature") %>%
  left_join(bulk_icgc %>% select(Signature, ICGC_Bulk = Std_Drift), by = "Signature") %>%
  mutate(
    # 拆解生理 80/20 癌旁中的纯肝细胞贡献 vs 非实质微环境增量
    Physio_Net_Crosstalk = Hep_80 - Hep_100
  ) %>%
  arrange(desc(Hep_80))

message("\n>>> [2/4] 生理比例校正 (80% Hep + 20% NonHep) vs 纯肝细胞 vs 真实 Bulk 对比表:")
print(as.data.frame(
  tri_wide %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(Signature,
           PureHep_100 = Hep_100,
           Physio_80_20 = Hep_80,
           NonHep_0 = Hep_0,
           TCGA_Bulk, ICGC_Bulk)
), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. 顺手核查 CYT_Rooney & T_cell_exhaustion 在纯肝细胞中正漂移的原因
# ------------------------------------------------------------------------------
check_genes <- intersect(c("GZMA", "PRF1", "TGFB1", "TGFBR1", "SMAD7", "SERPINE1", "ID1", "ID2", "ID3"), rownames(P_T))
chk_tbl <- tibble(
  Gene = check_genes,
  Tumor_PB_CPM   = round(rowMeans(P_T[check_genes, , drop = FALSE]), 1),
  PureHep_PB_CPM = round(rowMeans(P_Hep[check_genes, , drop = FALSE]), 1),
  NonHep_PB_CPM  = round(rowMeans(P_NonHep[check_genes, , drop = FALSE]), 1)
)
message("\n>>> [3/4] 典型基因在肿瘤 / 纯正常肝细胞 / 癌旁非实质中的平均 CPM:")
print(as.data.frame(chk_tbl), row.names = FALSE)

# 保存校正结果与升级版 Fig 3
saveRDS(
  list(sweep_df = sweep_df, sweep_cor = sweep_cor, tri_wide = tri_wide),
  file = "01_clean_data/03b_scrna_calibrated.rds"
)
write.csv(tri_wide, "03_results/03b_scrna_physiological_triarm.csv", row.names = FALSE)

# 绘制升级版 Figure 3 (左: 三臂生理校正对比 | 右: 肝细胞占比扫描与真实 Bulk 匹配曲线)
sig_ord <- rev(tri_wide$Signature)

p3a_cal <- tri_df %>%
  mutate(Signature = factor(Signature, levels = sig_ord)) %>%
  ggplot(aes(x = Std_Drift, y = Signature)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_line(aes(group = Signature), color = "grey72", linewidth = 0.85) +
  geom_point(aes(color = Arm_Label, shape = Arm_Label), size = 3.1) +
  scale_color_manual(
    values = c(
      "1_Pure_Hepatocytes (100% Hep)"                  = "#1F77B4",
      "2_Physiological_Liver (80% Hep + 20% NonHep)"   = "#D62728",
      "3_Pure_NonParenchymal (0% Hep, Immune/Stroma)"  = "#FF7F0E"
    ),
    name = "Admixture Composition"
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = "Admixture Composition") +
  labs(
    title    = "A. Physiologically Calibrated Single-Cell Mechanism Decomposition",
    subtitle = "30% admixture: Pure Hepatocytes vs. Calibrated Liver (80/20) vs. Pure Non-Parenchymal",
    x        = "Standardized Score Drift (ssGSEA, 30%)",
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        plot.title = element_text(face = "bold"))

p3b_cal <- sweep_cor %>%
  pivot_longer(cols = c(Rho_vs_TCGA, Rho_vs_ICGC),
               names_to = "Bulk_Cohort", values_to = "Spearman_Rho") %>%
  mutate(Bulk_Cohort = recode(Bulk_Cohort,
                              Rho_vs_TCGA = "vs. TCGA-LIHC Bulk",
                              Rho_vs_ICGC = "vs. ICGC-LIRI-JP Bulk")) %>%
  ggplot(aes(x = Hep_Frac, y = Spearman_Rho, color = Bulk_Cohort)) +
  annotate("rect", xmin = 70, xmax = 85, ymin = -Inf, ymax = Inf,
           fill = "#E8F5E9", alpha = 0.6) +
  annotate("text", x = 77.5, y = min(sweep_cor$Rho_vs_TCGA, sweep_cor$Rho_vs_ICGC) + 0.03,
           label = "Physiological\nliver mRNA\nrange (70-85%)", size = 3, color = "#2E7D32", fontface = "bold") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("vs. TCGA-LIHC Bulk" = "#B2182B", "vs. ICGC-LIRI-JP Bulk" = "#2166AC"),
                     name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  labs(
    title    = "B. Concordance with Real Bulk Admixture",
    subtitle = "Effect of correcting scRNA-seq hepatocyte dissociation bias",
    x        = "Hepatocyte RNA Fraction in Simulated Adjacent Liver (%)",
    y        = "Spearman rho (scRNA-seq vs. Bulk 30% Drift)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p3_cal <- (p3a_cal | p3b_cal) + plot_layout(widths = c(1.35, 1))
ggsave("04_figures/Fig3_Calibrated_scRNA_Mechanism.pdf", p3_cal, width = 13.5, height = 8)
ggsave("04_figures/Fig3_Calibrated_scRNA_Mechanism.png", p3_cal, width = 13.5, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 5. 预览 GSE238264 空间压缩包内部文件名（不解压，为第 4 关铺路）
# ------------------------------------------------------------------------------
f_sp1 <- "00_raw_data/spatial/GSE238264/GSE238264_RAW/GSM7661255_HCC1R.tar.gz"
message("\n>>> [4/4] 预览 GSE238264 单个样本压缩包内部结构:")
if (file.exists(f_sp1)) {
  print(head(untar(f_sp1, list = TRUE), 25))
} else {
  message("    未找到 ", f_sp1)
}
message(">>> 全部完成！请查看 04_figures/Fig3_Calibrated_scRNA_Mechanism.png。")
