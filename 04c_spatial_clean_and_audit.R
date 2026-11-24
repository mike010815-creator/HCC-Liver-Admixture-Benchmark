# ==============================================================================
# 脚本名称: 04c_spatial_clean_and_audit.R
# 核心任务:
#   1. 对空间切片 HCC3R 采用 AddModuleScore 扣除 UMI/测序深度混杂，重绘精修版 Fig 4
#   2. 跨模态闭环验证: 检验第 5 关构建的“单细胞净化特征 (S3)”在真实空间切片 HCC3R 上
#      是否同样显著降低了与 NHAS 的空间假性共变
#   3. 打印第 5 关 S3 净化特征在 789 例真实患者中的 OS 预后保真度对比表
#   4. 核验 Table 1 与全部 6 张主图的完整竣工清单
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Seurat)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 加载所需对象
# ------------------------------------------------------------------------------
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
bulk_sim <- readRDS("01_clean_data/02d_multicohort_sim_results.rds")
res4_cli <- readRDS("01_clean_data/04b_clinical_survival_results.rds")
res5_obj <- readRDS("01_clean_data/05_correction_benchmark_results.rds")

sig_list      <- sigs_obj$sig_list
sig_meta      <- sigs_obj$sig_meta
nhas_genes    <- res3_obj$nhas_genes
pur_sig_list  <- res5_obj$purified_sig_list
pruned_sigs   <- res5_obj$prune_summary %>% filter(Was_Pruned) %>% pull(Signature)

# ------------------------------------------------------------------------------
# 2. 读入强梯度代表切片 HCC3R，用 AddModuleScore 消除 UMI/坏死区深度混杂
# ------------------------------------------------------------------------------
h3_dir <- "01_clean_data/spatial_unpacked/GSE238264/HCC3R"
stopifnot("找不到已解压的 HCC3R 空间目录" = dir.exists(h3_dir))

message(">>> [1/4] 正在读入 HCC3R 空间切片并执行去 UMI 混杂 AddModuleScore 打分...")
cnt <- Read10X(file.path(h3_dir, "filtered_feature_bc_matrix"))
if (is.list(cnt)) cnt <- cnt[[1]]

seu <- CreateSeuratObject(counts = cnt, assay = "Spatial", project = "HCC3R")
img <- Read10X_Image(file.path(h3_dir, "spatial"),
                     image.name = "tissue_lowres_image.png",
                     filter.matrix = TRUE)
DefaultAssay(img) <- "Spatial"
common_bc <- intersect(colnames(seu), Cells(img))
seu <- subset(seu, cells = common_bc)
img <- img[common_bc]
seu[["slice1"]] <- img

# 过滤低质量坏死/空泡捕获点 (nFeature_Spatial >= 500)，进一步剔除技术空洞
seu <- subset(seu, subset = nFeature_Spatial >= 500)
seu <- NormalizeData(seu, verbose = FALSE)
message(sprintf("    质控后保留活组织 Spots: %d 个", ncol(seu)))

# 辅助函数: 安全执行 AddModuleScore 并返回每个 Spot 的去背景得分向量
score_module_safe <- function(obj, gene_vec, col_prefix) {
  gu <- intersect(gene_vec, rownames(obj))
  if (length(gu) < 2) return(rep(NA_real_, ncol(obj)))
  tmp <- AddModuleScore(obj, features = list(gu), name = col_prefix, ctrl = 50, seed = 2026)
  as.numeric(tmp[[paste0(col_prefix, "1")]][, 1])
}

seu$NHAS_MS <- score_module_safe(seu, nhas_genes, "TMP_NHAS_")

# 计算全部 24 个原始特征的 AddModuleScore
raw_ms_mat <- sapply(names(sig_list), function(s) {
  score_module_safe(seu, sig_list[[s]], "TMP_RAW_")
})

# 同时计算 18 个 S3 单细胞净化特征的 AddModuleScore，检验空间去混杂闭环
pur_ms_mat <- sapply(pruned_sigs, function(s) {
  score_module_safe(seu, pur_sig_list[[s]], "TMP_PUR_")
})

# 把典型特征写入 Seurat metadata 供空间热图绘制
for (ft in c("CTRL_G2M_Checkpoint", "HM_Complement", "T_cell_exhaustion",
             "HM_TGF_beta", "TAM_Peng", "MDSC_Peng")) {
  seu[[paste0("MS_", ft)]] <- raw_ms_mat[, ft]
}
seu$MS_T_cell_exhaustion_Purified <- pur_ms_mat[, "T_cell_exhaustion"]
seu$MS_HM_Complement_Purified     <- pur_ms_mat[, "HM_Complement"]

# ------------------------------------------------------------------------------
# 3. 量化去混杂后的空间相关性，以及 S3 净化特征在空间切片上的减负效果
# ------------------------------------------------------------------------------
bulk_mean30 <- bulk_sim %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  group_by(Signature) %>%
  summarise(Bulk30_Mean = mean(Std_Drift), .groups = "drop")

sp_clean_df <- tibble(
  Signature = colnames(raw_ms_mat),
  Spatial_MS_Rho_Raw = apply(raw_ms_mat, 2, function(v) {
    cor(seu$NHAS_MS, v, method = "spearman", use = "complete.obs")
  })
) %>%
  left_join(
    tibble(
      Signature = colnames(pur_ms_mat),
      Spatial_MS_Rho_Pur = apply(pur_ms_mat, 2, function(v) {
        cor(seu$NHAS_MS, v, method = "spearman", use = "complete.obs")
      })
    ),
    by = "Signature"
  ) %>%
  left_join(bulk_mean30, by = "Signature") %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature") %>%
  arrange(desc(Spatial_MS_Rho_Raw))

r_clean_bulk <- cor(sp_clean_df$Spatial_MS_Rho_Raw, sp_clean_df$Bulk30_Mean, method = "spearman")
message(sprintf("\n>>> [2/4] 去 UMI 混杂后 HCC3R 空间梯度 vs 三队列 Bulk 模拟一致性: Spearman rho = %.3f",
                r_clean_bulk))
print(as.data.frame(
  sp_clean_df %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    select(Signature, Spatial_MS_Rho_Raw, Spatial_MS_Rho_Pur, Bulk30_Mean)
), row.names = FALSE)

# 统计 S3 净化在空间切片上是否显著压低了与 NHAS 的假性空间共变
sp_pur_cmp <- sp_clean_df %>% filter(!is.na(Spatial_MS_Rho_Pur))
w_test_sp  <- wilcox.test(sp_pur_cmp$Spatial_MS_Rho_Raw, sp_pur_cmp$Spatial_MS_Rho_Pur, paired = TRUE)
message(sprintf("    [空间闭环验证] 18 个特征净化前后空间 NHAS 相关性中位数: %.3f -> %.3f (配对 Wilcoxon P = %.3g)",
                median(sp_pur_cmp$Spatial_MS_Rho_Raw),
                median(sp_pur_cmp$Spatial_MS_Rho_Pur),
                w_test_sp$p.value))

write.csv(sp_clean_df, "03_results/04c_hcc3r_spatial_deconfounded.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. 重绘精修版 Figure 4 (上: 6 张空间彩图 | 下左: 空间共变 vs Bulk | 下右: 空间净化效果)
# ------------------------------------------------------------------------------
message("\n>>> [3/4] 正在绘制去混杂精修版 Figure 4...")

# 上半部分: 6 张代表性空间分布图 (明确呈现正相关 vs 负相关增殖对照)
map_specs <- c(
  "NHAS_MS"               = "NHAS (Non-tumor Hepatocyte)",
  "MS_CTRL_G2M_Checkpoint" = "CTRL_G2M_Checkpoint (Neg. Ctrl)",
  "MS_HM_Complement"      = "HM_Complement (High Vuln.)",
  "MS_T_cell_exhaustion"  = "T_cell_exhaustion (Raw)",
  "MS_T_cell_exhaustion_Purified" = "T_cell_exhaustion (Purified)",
  "MS_TAM_Peng"           = "TAM_Peng (Myeloid)"
)

sp_plot_list <- lapply(names(map_specs), function(col_n) {
  SpatialFeaturePlot(seu, features = col_n, image.alpha = 0.35, pt.size.factor = 1.6) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = median(seu[[col_n]][, 1], na.rm = TRUE),
                         name = "Score") +
    ggtitle(map_specs[[col_n]]) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.key.height = unit(0.45, "cm"))
})
p4_top <- wrap_plots(sp_plot_list, ncol = 3)

# 下左: 24 个特征去混杂空间共变条形图 (一目了然: 补体/炎症/耗竭在右侧正值，增殖对照在左侧负值)
sig_ord4c <- rev(sp_clean_df$Signature)
p4_bl <- sp_clean_df %>%
  mutate(
    Signature = factor(Signature, levels = sig_ord4c),
    Dir = ifelse(Spatial_MS_Rho_Raw >= 0, "Positive (Co-enriched with NHAS)", "Negative (Tumor-Core Enriched)")
  ) %>%
  ggplot(aes(x = Spatial_MS_Rho_Raw, y = Signature, fill = Dir)) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
  geom_col(width = 0.72, alpha = 0.88) +
  scale_fill_manual(values = c("Positive (Co-enriched with NHAS)" = "#D62728",
                               "Negative (Tumor-Core Enriched)"    = "#1F77B4"),
                    name = NULL) +
  labs(
    title    = "B. Deconfounded Spatial Co-variation with NHAS (HCC3R)",
    subtitle = sprintf("Background-adjusted AddModuleScore (vs. Bulk 30%% Drift: rho = %.3f)", r_clean_bulk),
    x        = "Spatial Spearman rho with NHAS", y = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 下右: 策略 3 基因级净化在真实空间切片 HCC3R 上的去共变效果对比
pur_ord4c <- sp_pur_cmp %>% arrange(Spatial_MS_Rho_Raw) %>% pull(Signature)
p4_br <- sp_pur_cmp %>%
  select(Signature, Spatial_MS_Rho_Raw, Spatial_MS_Rho_Pur) %>%
  pivot_longer(-Signature, names_to = "Version", values_to = "Spatial_Rho") %>%
  mutate(
    Signature = factor(Signature, levels = pur_ord4c),
    Version   = recode(Version,
                       Spatial_MS_Rho_Raw = "Original Signature (Raw)",
                       Spatial_MS_Rho_Pur = "scRNA-Purified Signature (S3)")
  ) %>%
  ggplot(aes(x = Spatial_Rho, y = Signature)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_line(aes(group = Signature), color = "grey65", linewidth = 0.85) +
  geom_point(aes(color = Version), size = 2.8) +
  scale_color_manual(values = c("Original Signature (Raw)"      = "#D62728",
                                "scRNA-Purified Signature (S3)" = "#2CA02C"),
                     name = NULL) +
  labs(
    title    = "C. Spatial Deconfounding by scRNA Gene Purification",
    subtitle = sprintf("Paired Wilcoxon P = %.2g across 18 pruned signatures", w_test_sp$p.value),
    x        = "Spatial Spearman rho with NHAS (HCC3R)", y = NULL
  ) +
  theme_bw(base_size = 9.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p4_clean <- p4_top / (p4_bl | p4_br) + plot_layout(heights = c(1.05, 1.25))
ggsave("04_figures/Fig4_Spatial_Clean_Validation.pdf", p4_clean, width = 14.5, height = 13.5)
ggsave("04_figures/Fig4_Spatial_Clean_Validation.png", p4_clean, width = 14.5, height = 13.5, dpi = 300)

# ------------------------------------------------------------------------------
# 5. 打印策略 3 净化特征的 789 例真实总生存预后保真度对比 + 全项目文件验收
# ------------------------------------------------------------------------------
message("\n>>> [4/4] 策略 3 净化特征在 789 例真实患者中的总生存 (OS) 预后保真度与 NHAS 去混杂核查:")
os_audit <- res4_cli$meta_cox_df %>%
  select(Signature, RealRho_Raw = Mean_Real_Rho, HR_Raw, P_Raw) %>%
  inner_join(res5_obj$pur_meta_cox, by = "Signature") %>%
  filter(Signature %in% pruned_sigs) %>%
  mutate(across(c(RealRho_Raw, Mean_RealRho_Pur, HR_Raw, HR_Purified_Meta), ~ round(.x, 3)),
         across(c(P_Raw, P_Purified_Meta), ~ signif(.x, 2))) %>%
  arrange(P_Raw)

print(as.data.frame(os_audit), row.names = FALSE)

message("\n>>> [全套论文图表与总表竣工清单]:")
print(list.files("04_figures", pattern = "\\.png$"))
print(list.files("03_results", pattern = "Table1"))
message(">>> 全部完成！请查看 04_figures/Fig4_Spatial_Clean_Validation.png 并把控制台输出发给我。")
