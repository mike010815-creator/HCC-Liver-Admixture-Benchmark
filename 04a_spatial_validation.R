# ==============================================================================
# 脚本名称: 04a_spatial_validation.R
# 核心任务:
#   1. 自动解压 GSE238264 全部 7 例 HCC 10x Visium 空间转录组切片
#   2. 逐张切片计算每个 Spot 的 NHAS (22 基因) 与 24 个免疫特征得分 (零 hdf5r 依赖)
#   3. 在每张切片上按 NHAS 划分空间区域 (Core / Transition / High-NHAS Margin)
#      并跑区域 Pseudo-bulk ssGSEA，检验真实空间梯度与 Bulk 模拟敏感性的一致性
#   4. 导出代表性切片空间热图 + 空间梯度验证总图 (Fig4_Spatial_Validation)
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Seurat)
  library(GSVA)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 加载免疫特征库、22 基因 NHAS 清单与 Bulk 模拟结果
# ------------------------------------------------------------------------------
sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
res3_obj <- readRDS("01_clean_data/03_scrna_mechanism_and_nhas.rds")
bulk_sim <- readRDS("01_clean_data/02d_multicohort_sim_results.rds")

sig_list   <- sigs_obj$sig_list
sig_meta   <- sigs_obj$sig_meta
nhas_genes <- res3_obj$nhas_genes

# ------------------------------------------------------------------------------
# 2. 自动解压 GSE238264 的 7 个 .tar.gz 样本包
# ------------------------------------------------------------------------------
raw_sp_dir  <- "00_raw_data/spatial/GSE238264/GSE238264_RAW"
unpack_root <- "01_clean_data/spatial_unpacked/GSE238264"
dir.create(unpack_root, recursive = TRUE, showWarnings = FALSE)

tar_files <- sort(list.files(raw_sp_dir, pattern = "\\.tar\\.gz$", full.names = TRUE))
stopifnot("未在 GSE238264_RAW 中找到 .tar.gz 文件" = length(tar_files) > 0)

message(sprintf(">>> [1/5] 检测到 %d 个 Visium 压缩包，检查/执行解压...", length(tar_files)))
for (tf in tar_files) {
  # 提取样本文件夹名（如 HCC1R）
  top_dir <- sub("/.*$", "", untar(tf, list = TRUE)[1])
  if (!dir.exists(file.path(unpack_root, top_dir, "spatial"))) {
    untar(tf, exdir = unpack_root)
  }
}

# 列出所有包含 spatial/ 子目录的有效样本文件夹
cand_dirs <- list.dirs(unpack_root, recursive = FALSE, full.names = TRUE)
sample_dirs <- cand_dirs[file.exists(file.path(cand_dirs, "spatial", "scalefactors_json.json"))]
sample_names <- basename(sample_dirs)
message("    已就绪切片列表: ", paste(sample_names, collapse = ", "))

# ------------------------------------------------------------------------------
# 3. 逐张读取切片（兼容 Seurat v4/v5，零 hdf5r 依赖），计算 Spot 级与区域级指标
# ------------------------------------------------------------------------------
get_layer_safe <- function(seu, layer_name) {
  tryCatch(
    GetAssayData(seu, assay = "Spatial", layer = layer_name),
    error = function(e) GetAssayData(seu, assay = "Spatial", slot = layer_name)
  )
}

slide_qc_list      <- list()
spot_cor_list      <- list()
zone_ssgsea_list   <- list()
seu_objects_cache  <- list()

message("\n>>> [2/5] 正在逐张处理 7 例 Visium 切片并量化空间 NHAS 梯度...")

for (i in seq_along(sample_dirs)) {
  s_dir  <- sample_dirs[i]
  s_name <- sample_names[i]
  
  # 读取计数矩阵（优先三件套文件夹，无需 hdf5r）
  mtx_dir <- file.path(s_dir, "filtered_feature_bc_matrix")
  counts  <- if (dir.exists(mtx_dir)) {
    Read10X(mtx_dir)
  } else {
    Read10X_h5(file.path(s_dir, "filtered_feature_bc_matrix.h5"))
  }
  if (is.list(counts)) counts <- counts[[1]]
  
  # 构建 Seurat 对象并挂载空间图像与坐标
  seu <- CreateSeuratObject(counts = counts, assay = "Spatial", project = s_name)
  img <- Read10X_Image(
    image.dir     = file.path(s_dir, "spatial"),
    image.name    = "tissue_lowres_image.png",
    filter.matrix = TRUE
  )
  DefaultAssay(img) <- "Spatial"
  common_bc <- intersect(colnames(seu), Cells(img))
  seu <- subset(seu, cells = common_bc)
  img <- img[common_bc]
  seu[["slice1"]] <- img
  
  # 基础质控过滤极低计数空点，再做 LogNormalize
  seu <- subset(seu, subset = nFeature_Spatial >= 200)
  seu <- NormalizeData(seu, assay = "Spatial", normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  
  cnt_mat  <- get_layer_safe(seu, "counts")
  norm_mat <- get_layer_safe(seu, "data")
  
  # Spot 级 NHAS 得分 (22 个成熟正常肝细胞基因的平均 LogNormalized 表达)
  nhas_use <- intersect(nhas_genes, rownames(norm_mat))
  nhas_spot <- Matrix::colMeans(norm_mat[nhas_use, , drop = FALSE])
  seu$NHAS_Spot <- nhas_spot
  
  # 同时计算 24 个免疫/对照特征的 Spot 级得分 (平均 LogNormalized 表达)
  sig_spot_mat <- sapply(names(sig_list), function(sn) {
    gu <- intersect(sig_list[[sn]], rownames(norm_mat))
    if (length(gu) < 2) return(rep(NA_real_, ncol(norm_mat)))
    Matrix::colMeans(norm_mat[gu, , drop = FALSE])
  })
  # 将关键展示特征写入 Seurat metadata 方便后续画空间热图
  for (key_f in c("T_cell_exhaustion", "HM_Complement", "HM_TGF_beta",
                  "HM_TNFA_via_NFKB", "T_cell_inflamed_GEP", "MDSC_Peng", "CTRL_G2M_Checkpoint")) {
    if (key_f %in% colnames(sig_spot_mat)) {
      seu[[key_f]] <- sig_spot_mat[, key_f]
    }
  }
  
  # 切片 NHAS 空间异质性评估 (识别哪些切片同时包含清晰肿瘤区与高 NHAS 癌旁/交界区)
  q_nhas <- quantile(nhas_spot, probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
  slide_qc_list[[s_name]] <- tibble(
    Slide           = s_name,
    N_Spots         = ncol(seu),
    NHAS_Detected   = length(nhas_use),
    NHAS_P10        = round(q_nhas[1], 3),
    NHAS_Median     = round(q_nhas[3], 3),
    NHAS_P90        = round(q_nhas[5], 3),
    NHAS_Range_90_10 = round(q_nhas[5] - q_nhas[1], 3),
    NHAS_SD         = round(sd(nhas_spot), 3),
    Pct_High_NHAS   = round(mean(nhas_spot > 1.0) * 100, 1)
  )
  
  # 指标 A: Spot 级连续空间相关性 Spearman(NHAS_Spot, Signature_Spot)
  spot_cor_list[[s_name]] <- tibble(
    Slide     = s_name,
    Signature = colnames(sig_spot_mat),
    Spot_Rho_with_NHAS = apply(sig_spot_mat, 2, function(v) {
      suppressWarnings(cor(nhas_spot, v, method = "spearman", use = "complete.obs"))
    })
  )
  
  # 指标 B: 将切片按 NHAS 四分位划分为 3 个空间带，聚合 Pseudo-bulk 跑严格 ssGSEA
  #   - Low_NHAS (Bottom 25%):  肿瘤核心区 (Tumor Core)
  #   - Mid_NHAS (25%–75%):     交界/过渡带 (Transition/Margin)
  #   - High_NHAS (Top 25%):    高混入/癌旁肝带 (High-Admixture / Adjacent)
  seu$NHAS_Zone <- case_when(
    nhas_spot <= q_nhas[2] ~ "1_Low_NHAS_Core",
    nhas_spot >= q_nhas[4] ~ "3_High_NHAS_Margin_Adj",
    TRUE                   ~ "2_Mid_NHAS_Transition"
  )
  
  zones <- c("1_Low_NHAS_Core", "2_Mid_NHAS_Transition", "3_High_NHAS_Margin_Adj")
  zone_cpm <- sapply(zones, function(zn) {
    idx <- which(seu$NHAS_Zone == zn)
    v <- Matrix::rowSums(cnt_mat[, idx, drop = FALSE])
    v / sum(v) * 1e6
  })
  zone_log2 <- log2(zone_cpm + 1)
  ss_param  <- ssgseaParam(as.matrix(zone_log2), sig_list, minSize = 2, normalize = FALSE)
  zone_ss   <- gsva(ss_param, verbose = FALSE)
  
  zone_ssgsea_list[[s_name]] <- as_tibble(zone_ss, rownames = "Signature") %>%
    mutate(
      Slide = s_name,
      Delta_Mid_minus_Low  = `2_Mid_NHAS_Transition`  - `1_Low_NHAS_Core`,
      Delta_High_minus_Low = `3_High_NHAS_Margin_Adj` - `1_Low_NHAS_Core`
    )
  
  seu_objects_cache[[s_name]] <- seu
  message(sprintf("    [%d/%d] %s 完成: %d spots | NHAS P10=%.2f -> P90=%.2f (差值 %.2f)",
                  i, length(sample_dirs), s_name, ncol(seu),
                  q_nhas[1], q_nhas[5], q_nhas[5] - q_nhas[1]))
}

slide_qc_df    <- bind_rows(slide_qc_list) %>% arrange(desc(NHAS_Range_90_10))
spot_cor_df    <- bind_rows(spot_cor_list)
zone_ssgsea_df <- bind_rows(zone_ssgsea_list)

message("\n>>> [3/5] 7 张 Visium 切片的 NHAS 空间梯度跨度质控表（按梯度从强到弱排序）:")
print(as.data.frame(slide_qc_df), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. 检验真实空间切片梯度 vs Bulk (TCGA / ICGC) 模拟漂移的跨模态一致性
# ------------------------------------------------------------------------------
# 取梯度跨度明显的切片（或全部 7 张切片取中位数）汇总每个特征的空间漂移方向
# 优先统计所有 7 张切片中位数 + Top 梯度切片（NHAS_Range_90_10 >= 0.5）中位数
strong_slides <- slide_qc_df %>% filter(NHAS_Range_90_10 >= 0.5) %>% pull(Slide)
if (length(strong_slides) < 3) strong_slides <- head(slide_qc_df$Slide, 4)

spatial_summary <- spot_cor_df %>%
  group_by(Signature) %>%
  summarise(
    All7_Median_Spot_Rho   = median(Spot_Rho_with_NHAS, na.rm = TRUE),
    Strong_Median_Spot_Rho = median(Spot_Rho_with_NHAS[Slide %in% strong_slides], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    zone_ssgsea_df %>%
      group_by(Signature) %>%
      summarise(
        Median_Zone_ssGSEA_Delta = median(Delta_High_minus_Low, na.rm = TRUE),
        Strong_Zone_ssGSEA_Delta = median(Delta_High_minus_Low[Slide %in% strong_slides], na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Signature"
  ) %>%
  left_join(
    bulk_sim %>%
      filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
      mutate(C_short = ifelse(grepl("TCGA", Cohort), "TCGA",
                              ifelse(grepl("ICGC", Cohort), "ICGC", "GSE14520"))) %>%
      select(Signature, C_short, Std_Drift) %>%
      pivot_wider(names_from = C_short, values_from = Std_Drift, names_prefix = "Bulk30_"),
    by = "Signature"
  ) %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature") %>%
  arrange(desc(Strong_Median_Spot_Rho))

r_sp_tcga_spot <- cor(spatial_summary$Strong_Median_Spot_Rho, spatial_summary$Bulk30_TCGA, method = "spearman")
r_sp_icgc_spot <- cor(spatial_summary$Strong_Median_Spot_Rho, spatial_summary$Bulk30_ICGC, method = "spearman")
r_sp_tcga_zone <- cor(spatial_summary$Strong_Zone_ssGSEA_Delta, spatial_summary$Bulk30_TCGA, method = "spearman")
r_sp_icgc_zone <- cor(spatial_summary$Strong_Zone_ssGSEA_Delta, spatial_summary$Bulk30_ICGC, method = "spearman")

message(sprintf("\n>>> [4/5] 真实空间切片 (n=%d 张强梯度切片: %s) vs Bulk 模拟跨模态一致性:",
                length(strong_slides), paste(strong_slides, collapse = ", ")))
message(sprintf("    [Spot 级 Spearman(NHAS, Sig)] vs TCGA Bulk 30%%: rho = %.3f | vs ICGC Bulk 30%%: rho = %.3f",
                r_sp_tcga_spot, r_sp_icgc_spot))
message(sprintf("    [空间区域 Pseudo-bulk ssGSEA]  vs TCGA Bulk 30%%: rho = %.3f | vs ICGC Bulk 30%%: rho = %.3f",
                r_sp_tcga_zone, r_sp_icgc_zone))

message("\n>>> 24 个特征真实空间梯度表现 vs Bulk 模拟对照表:")
print(as.data.frame(
  spatial_summary %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    select(Signature, Strong_Spot_Rho = Strong_Median_Spot_Rho,
           Strong_Zone_Delta = Strong_Zone_ssGSEA_Delta,
           Bulk30_TCGA, Bulk30_ICGC)
), row.names = FALSE)

# 保存空间统计结果（不把庞大的 Seurat 对象全存进单个大文件，只存轻量统计结果与代表切片）
best_slide_name <- slide_qc_df$Slide[1]
saveRDS(
  list(
    slide_qc_df     = slide_qc_df,
    spot_cor_df     = spot_cor_df,
    zone_ssgsea_df  = zone_ssgsea_df,
    spatial_summary = spatial_summary,
    strong_slides   = strong_slides,
    best_slide_name = best_slide_name
  ),
  file = "01_clean_data/04a_spatial_validation_results.rds"
)
write.csv(slide_qc_df,     "03_results/04a_spatial_slide_nhas_qc.csv", row.names = FALSE)
write.csv(spatial_summary, "03_results/04a_spatial_vs_bulk_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. 绘制 Figure 4: 空间热图展示 + 全切片梯度热图 + 空间 vs Bulk 跨模态验证散点图
# ------------------------------------------------------------------------------
message(sprintf("\n>>> [5/5] 正在绘制 Figure 4（代表切片: %s）...", best_slide_name))

# 图 4A: 全部 7 张切片的 Spot 级 NHAS 相关性热图
sig_ord4 <- rev(spatial_summary$Signature)
p4_heat <- spot_cor_df %>%
  mutate(
    Signature = factor(Signature, levels = sig_ord4),
    Slide     = factor(Slide, levels = slide_qc_df$Slide)
  ) %>%
  ggplot(aes(x = Slide, y = Signature, fill = Spot_Rho_with_NHAS)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Spot_Rho_with_NHAS)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "Spot rho\nwith NHAS") +
  labs(
    title    = "A. Spatial Co-variation with Non-tumor Hepatocyte Score (NHAS)",
    subtitle = "7 HCC Visium slides (GSE238264), ordered left-to-right by NHAS gradient width",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
        plot.title = element_text(face = "bold"))

# 图 4B: 空间区域 Pseudo-bulk ssGSEA 漂移 vs Bulk (ICGC) 30% 模拟漂移散点图
p4_sc <- ggplot(spatial_summary, aes(x = Bulk30_ICGC, y = Strong_Zone_ssGSEA_Delta, color = Axis)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.75, fill = "grey88") +
  geom_point(size = 2.8) +
  annotate("text", x = min(spatial_summary$Bulk30_ICGC), y = max(spatial_summary$Strong_Zone_ssGSEA_Delta),
           hjust = 0, vjust = 1, fontface = "bold", size = 3.5,
           label = sprintf("Spatial vs. ICGC Bulk: rho = %.3f\nSpatial vs. TCGA Bulk: rho = %.3f",
                           r_sp_icgc_zone, r_sp_tcga_zone)) +
  labs(
    title    = "B. Real Spatial Margin Drift vs. Bulk Simulation",
    subtitle = "High-NHAS vs. Low-NHAS spatial zones (ssGSEA)",
    x        = "In Silico Bulk 30% Admixture Drift (ICGC-LIRI-JP)",
    y        = "Real Spatial Zone Drift (High-NHAS minus Low-NHAS)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        plot.title = element_text(face = "bold"))

# 图 4C: 代表性切片（NHAS 跨度最大）三区梯度箱线图（精选 6 个对比特征）
best_seu <- seu_objects_cache[[best_slide_name]]
key_show <- c("T_cell_exhaustion", "HM_Complement", "HM_TGF_beta",
              "T_cell_inflamed_GEP", "MDSC_Peng", "CTRL_G2M_Checkpoint")
key_show <- intersect(key_show, colnames(best_seu@meta.data))

box_df <- best_seu@meta.data %>%
  select(NHAS_Zone, all_of(key_show)) %>%
  pivot_longer(-NHAS_Zone, names_to = "Signature", values_to = "Spot_Score") %>%
  mutate(
    Signature = factor(Signature, levels = key_show),
    Zone_Short = recode(NHAS_Zone,
                        "1_Low_NHAS_Core"        = "Core\n(Low NHAS)",
                        "2_Mid_NHAS_Transition"  = "Transition\n(Mid NHAS)",
                        "3_High_NHAS_Margin_Adj" = "Margin/Adj\n(High NHAS)")
  )

p4_box <- ggplot(box_df, aes(x = Zone_Short, y = Spot_Score, fill = Zone_Short)) +
  geom_violin(alpha = 0.45, color = NA) +
  geom_boxplot(width = 0.25, outlier.size = 0.4, alpha = 0.85) +
  facet_wrap(~ Signature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("#3B4992", "#E69F00", "#D62728"), guide = "none") +
  labs(
    title    = sprintf("C. Representative Slide (%s): Gradient Across Spatial NHAS Zones", best_slide_name),
    x        = NULL, y = "Spot-Level Signature Score"
  ) +
  theme_bw(base_size = 9.5) +
  theme(strip.background = element_rect(fill = "#F2F4F7"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

p4_main <- (p4_heat | (p4_sc / p4_box)) + plot_layout(widths = c(1.05, 1.35))
ggsave("04_figures/Fig4_Spatial_Gradient_Validation.pdf", p4_main, width = 15, height = 9)
ggsave("04_figures/Fig4_Spatial_Gradient_Validation.png", p4_main, width = 15, height = 9, dpi = 300)

# 额外导出代表性切片的真实组织空间分布彩图 (Fig4S_Representative_Spatial_Maps.png)
tryCatch({
  sp_feats <- c("NHAS_Spot", "CTRL_G2M_Checkpoint", "T_cell_exhaustion",
                "HM_Complement", "HM_TGF_beta", "MDSC_Peng")
  sp_feats <- intersect(sp_feats, colnames(best_seu@meta.data))
  sp_plots <- lapply(sp_feats, function(ft) {
    SpatialFeaturePlot(best_seu, features = ft, image.alpha = 0.45) +
      ggtitle(ft) +
      theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
            legend.position = "right")
  })
  p4_maps <- wrap_plots(sp_plots, ncol = 3)
  ggsave("04_figures/Fig4S_Representative_Spatial_Maps.pdf", p4_maps, width = 13.5, height = 8.5)
  ggsave("04_figures/Fig4S_Representative_Spatial_Maps.png", p4_maps, width = 13.5, height = 8.5, dpi = 300)
}, error = function(e) {
  message("    提示: SpatialFeaturePlot 渲染遇到小异常（不影响主图 Fig4）: ", e$message)
})

message(">>> 全部完成！请查看 04_figures/Fig4_Spatial_Gradient_Validation.png 与 Fig4S_Representative_Spatial_Maps.png。")
