# ==============================================================================
# 脚本名称: 03_scrna_pseudobulk.R
# 核心任务:
#   1. 高效分块读取 GSE149614 单细胞计数矩阵并转为稀疏矩阵缓存 (.rds)
#   2. 在配对原发癌 (Tumor) + 癌旁肝 (Normal) 患者中构造双臂 Pseudo-bulk 模拟:
#      - Arm A: 纯正常肝细胞混入 (Hepatocytes only)
#      - Arm B: 完整癌旁非肿瘤肝混入 (Full Adjacent Liver)
#   3. 机制拆解 + 锁定每个高敏感特征的“污染易感基因”与主导细胞来源 (产出 Fig 3)
#   4. 构建零循环论证的正常肝实质混入评分 (NHAS, 30 基因) 并在三队列验证精度
# ==============================================================================

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Matrix)
  library(GSVA)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# 1. 读取并缓存 GSE149614 稀疏计数矩阵与细胞注释表
# ------------------------------------------------------------------------------
cache_file <- "01_clean_data/gse149614_sparse.rds"

if (file.exists(cache_file)) {
  message(">>> [1/6] 检测到已缓存的稀疏单细胞对象，秒速加载中...")
  sc_obj <- readRDS(cache_file)
  counts_sp <- sc_obj$counts
  meta      <- sc_obj$meta
} else {
  f_meta <- "00_raw_data/scRNA_GSE149614/GSE149614_HCC.metadata.updated.txt.gz"
  f_cnt  <- "00_raw_data/scRNA_GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
  stopifnot(file.exists(f_meta), file.exists(f_cnt))
  
  message(">>> [1/6] 首次读取 GSE149614（约 7.2 万细胞，耗时 1–2 分钟，请耐心等待）...")
  meta_raw <- fread(f_meta, data.table = FALSE)
  cnt_dt   <- fread(f_cnt, data.table = TRUE)
  
  # 第一列为基因名，其余列为细胞 ID
  gene_names <- as.character(cnt_dt[[1]])
  cell_cols  <- colnames(cnt_dt)[-1]
  
  # 自动识别 meta 中与计数矩阵列名匹配的细胞 ID 列
  overlap_counts <- sapply(meta_raw, function(col) sum(as.character(col) %in% cell_cols))
  id_col_idx <- which.max(overlap_counts)
  stopifnot("metadata 中未匹配到细胞 ID 列" = overlap_counts[id_col_idx] > 1000)
  
  meta_raw$cell_id <- as.character(meta_raw[[id_col_idx]])
  common_cells <- intersect(cell_cols, meta_raw$cell_id)
  message(sprintf("    成功匹配细胞数: %d | 基因数: %d", length(common_cells), length(gene_names)))
  
  # 仅保留匹配到的细胞，并按 5000 细胞分块转为稀疏矩阵 dgCMatrix（节省内存、防溢出）
  meta <- meta_raw[match(common_cells, meta_raw$cell_id), , drop = FALSE]
  rownames(meta) <- meta$cell_id
  
  chunk_size <- 5000L
  starts <- seq(1L, length(common_cells), by = chunk_size)
  sp_list <- vector("list", length(starts))
  
  for (i in seq_along(starts)) {
    e <- min(starts[i] + chunk_size - 1L, length(common_cells))
    sub_cells <- common_cells[starts[i]:e]
    sub_mat <- as.matrix(cnt_dt[, ..sub_cells])
    storage.mode(sub_mat) <- "double"
    sp_list[[i]] <- Matrix(sub_mat, sparse = TRUE)
    if (i %% 5 == 0 || i == length(starts)) {
      message(sprintf("    稀疏化进度: %d / %d 细胞", e, length(common_cells)))
    }
  }
  counts_sp <- do.call(cbind, sp_list)
  rownames(counts_sp) <- gene_names
  
  # 去除重复基因名（保留总 UMI 最高的）
  if (any(duplicated(rownames(counts_sp)))) {
    ord_g <- order(rownames(counts_sp), -Matrix::rowSums(counts_sp))
    counts_sp <- counts_sp[ord_g, ]
    counts_sp <- counts_sp[!duplicated(rownames(counts_sp)), ]
  }
  
  rm(cnt_dt, sp_list); gc()
  saveRDS(list(counts = counts_sp, meta = meta), file = cache_file)
  message("    已缓存稀疏单细胞对象至 01_clean_data/gse149614_sparse.rds。")
}

# ------------------------------------------------------------------------------
# 2. 自动解析 patient / site / celltype 列并核对配对患者与细胞组成
# ------------------------------------------------------------------------------
message("\n>>> [2/6] 单细胞 Metadata 列名: ", paste(colnames(meta), collapse = ", "))

col_pat  <- grep("^patient", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
col_site <- grep("^site$|^tissue$|^origin$", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
col_ct   <- grep("^celltype$|^cell_type$|^major", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(col_pat), !is.na(col_site), !is.na(col_ct))

meta$Patient  <- as.character(meta[[col_pat]])
meta$Site     <- as.character(meta[[col_site]])
meta$CellType <- as.character(meta[[col_ct]])

message("    组织部位 x 细胞大类分布表:")
print(table(meta$Site, meta$CellType))

# 识别 Primary Tumor 与 Adjacent/Normal Liver 标签
is_T <- grepl("^Tumor$|^PT$|^Primary", meta$Site, ignore.case = TRUE)
is_N <- grepl("^Normal$|^NTL$|^Adj|^Non", meta$Site, ignore.case = TRUE)
is_Hep <- grepl("Hep", meta$CellType, ignore.case = TRUE)

# 筛选同时具有原发癌 (>=200 细胞) 和癌旁肝 (且含 >=50 个正常肝细胞) 的配对患者
pat_qc <- meta %>%
  group_by(Patient) %>%
  summarise(
    n_T       = sum(is_T[cur_group_rows()]),
    n_N_all   = sum(is_N[cur_group_rows()]),
    n_N_hep   = sum(is_N[cur_group_rows()] & is_Hep[cur_group_rows()]),
    n_N_nonhp = sum(is_N[cur_group_rows()] & !is_Hep[cur_group_rows()]),
    .groups   = "drop"
  )
message("\n    各患者细胞量质控表:")
print(as.data.frame(pat_qc), row.names = FALSE)

paired_pts <- pat_qc %>%
  filter(n_T >= 200, n_N_all >= 200, n_N_hep >= 50) %>%
  pull(Patient) %>% sort()

stopifnot("符合条件的配对单细胞患者不足 3 例" = length(paired_pts) >= 3)
message(sprintf("    入选双臂 Pseudo-bulk 配对患者共 %d 例: %s",
                length(paired_pts), paste(paired_pts, collapse = ", ")))

# ------------------------------------------------------------------------------
# 3. 构建每位患者的 3 套纯净 Pseudo-bulk 表达谱 (归一化为 1e6 CPM)
# ------------------------------------------------------------------------------
to_cpm <- function(mat_sp, cell_ids) {
  v <- Matrix::rowSums(mat_sp[, cell_ids, drop = FALSE])
  v / sum(v) * 1e6
}

P_T    <- sapply(paired_pts, function(p) to_cpm(counts_sp, meta$cell_id[meta$Patient == p & is_T]))
P_Hep  <- sapply(paired_pts, function(p) to_cpm(counts_sp, meta$cell_id[meta$Patient == p & is_N & is_Hep]))
P_NTL  <- sapply(paired_pts, function(p) to_cpm(counts_sp, meta$cell_id[meta$Patient == p & is_N]))

# 同时构造纯“癌旁非肝细胞”谱 (用于校验纯微环境贡献)
P_NTL_nonhep <- sapply(paired_pts, function(p) {
  to_cpm(counts_sp, meta$cell_id[meta$Patient == p & is_N & !is_Hep])
})

sigs_obj <- readRDS("01_clean_data/immune_signatures.rds")
sig_list <- sigs_obj$sig_list
sig_meta <- sigs_obj$sig_meta

# 运行双臂 0%–50% 混合模拟 (ssGSEA 与 Zscore)
alphas <- c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

run_pb_arm <- function(N_ref, arm_name) {
  mix_logs <- lapply(alphas, function(a) log2((1 - a) * P_T + a * N_ref + 1))
  
  # ssGSEA
  ss_list <- lapply(mix_logs, function(m) {
    p <- ssgseaParam(as.matrix(m), sig_list, minSize = 2, normalize = FALSE)
    gsva(p, verbose = FALSE)
  })
  # Fixed-baseline Zscore
  m0 <- mix_logs[[1]]
  mu0 <- rowMeans(m0); sd0 <- apply(m0, 1, sd); sd0[sd0 == 0] <- NA
  zs_list <- lapply(mix_logs, function(m) {
    zm <- sweep(sweep(m, 1, mu0, "-"), 1, sd0, "/")
    out <- sapply(sig_list, function(g) {
      gu <- intersect(g, rownames(zm))
      colMeans(zm[gu, , drop = FALSE], na.rm = TRUE)
    })
    t(out)
  })
  
  calc_arm <- function(slist, method_lbl) {
    s0_mat <- slist[[1]]
    map_dfr(seq_along(alphas), function(k) {
      a <- alphas[k]; sa_mat <- slist[[k]]
      map_dfr(rownames(s0_mat), function(sig) {
        s0 <- as.numeric(s0_mat[sig, ]); sa <- as.numeric(sa_mat[sig, ])
        sd0_val <- sd(s0)
        diffs <- (sa - s0) / sd0_val
        tibble(
          Arm       = arm_name,
          Method    = method_lbl,
          Signature = sig,
          Alpha_Pct = a * 100,
          Std_Drift = mean(diffs),
          SEM_Drift = sd(diffs) / sqrt(length(diffs))
        )
      })
    })
  }
  bind_rows(calc_arm(ss_list, "ssGSEA"), calc_arm(zs_list, "Zscore"))
}

message("\n>>> [3/6] 正在运行双臂单细胞 Pseudo-bulk 控制变量混合模拟...")
pb_res <- bind_rows(
  run_pb_arm(P_Hep, "ArmA_Pure_Hepatocytes"),
  run_pb_arm(P_NTL, "ArmB_Full_Adjacent_Liver")
) %>%
  left_join(sig_meta %>% select(Signature, Axis), by = "Signature")

# 30% 混入下的机制拆解表 (Arm A = 肝细胞效应, Arm B - Arm A = 癌旁非实质细胞串扰效应)
decomp_30 <- pb_res %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  select(Signature, Axis, Arm, Std_Drift) %>%
  pivot_wider(names_from = Arm, values_from = Std_Drift) %>%
  mutate(
    Hepatocyte_Component    = ArmA_Pure_Hepatocytes,
    NonParenchymal_Crosstalk = ArmB_Full_Adjacent_Liver - ArmA_Pure_Hepatocytes,
    Total_NTL_Drift         = ArmB_Full_Adjacent_Liver,
    Dominant_Mechanism = case_when(
      Total_NTL_Drift < -0.20 ~ "Type_III_Parenchymal_Dilution",
      Hepatocyte_Component >= 0.25 & Hepatocyte_Component >= NonParenchymal_Crosstalk ~
        "Type_I_Hepatocyte_Intrinsic_Intrusion",
      NonParenchymal_Crosstalk >= 0.25 & NonParenchymal_Crosstalk > Hepatocyte_Component ~
        "Type_II_Adjacent_Microenv_Crosstalk",
      abs(Total_NTL_Drift) < 0.25 ~ "Relatively_Resistant",
      TRUE ~ "Mixed_Hep_and_Microenv"
    )
  ) %>%
  arrange(desc(Total_NTL_Drift))

message("\n>>> [机制拆解表] 30% 混入 (ssGSEA): 纯肝细胞效应 (Arm A) vs 完整癌旁效应 (Arm B):")
print(as.data.frame(
  decomp_30 %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(Signature,
           ArmA_Hep = Hepatocyte_Component,
           ArmB_FullNTL = Total_NTL_Drift,
           Crosstalk_BminusA = NonParenchymal_Crosstalk,
           Dominant_Mechanism)
), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. 基因级“罪魁祸首”溯源（精确到每个癌旁细胞亚群的平均 CPM 与检出率）
# ------------------------------------------------------------------------------
message("\n>>> [4/6] 正在执行基因级细胞来源溯源（定位污染易感基因）...")

# 仅保留 Primary Tumor 与 Adjacent Normal 的细胞
tn_idx <- which(is_T | is_N)
sp_tn  <- counts_sp[, tn_idx]
meta_tn <- meta[tn_idx, ]
meta_tn$Compartment <- ifelse(
  grepl("^Tumor$|^PT$|^Primary", meta_tn$Site, ignore.case = TRUE),
  paste0("T_", meta_tn$CellType),
  paste0("N_", meta_tn$CellType)
)

# 对单细胞做每细胞 CPM 归一化（稀疏矩阵列归一化，几秒完成）
lib_sizes <- Matrix::colSums(sp_tn)
sp_tn_cpm <- sp_tn %*% Diagonal(x = 1e6 / pmax(lib_sizes, 1))

compartments <- sort(unique(meta_tn$Compartment))
mean_cpm_by_comp <- sapply(compartments, function(cp) {
  idx <- which(meta_tn$Compartment == cp)
  Matrix::rowMeans(sp_tn_cpm[, idx, drop = FALSE])
})
pct_by_comp <- sapply(compartments, function(cp) {
  idx <- which(meta_tn$Compartment == cp)
  Matrix::rowMeans(sp_tn[, idx, drop = FALSE] > 0)
})

# 同时计算癌旁全组织 vs 原发癌全组织的单细胞伪整体均值比值 (log2FC)
all_sig_genes <- intersect(unique(unlist(sig_list)), rownames(sp_tn))
nt_comp_cols  <- grep("^N_", colnames(mean_cpm_by_comp), value = TRUE)
t_comp_cols   <- grep("^T_", colnames(mean_cpm_by_comp), value = TRUE)

gene_source_df <- tibble(
  Gene               = all_sig_genes,
  PB_Tumor_CPM       = rowMeans(P_T[all_sig_genes, , drop = FALSE]),
  PB_NormHep_CPM     = rowMeans(P_Hep[all_sig_genes, , drop = FALSE]),
  PB_FullNTL_CPM     = rowMeans(P_NTL[all_sig_genes, , drop = FALSE]),
  Log2FC_NTL_vs_T    = log2((PB_FullNTL_CPM + 1) / (PB_Tumor_CPM + 1)),
  Log2FC_NormHep_vs_T = log2((PB_NormHep_CPM + 1) / (PB_Tumor_CPM + 1)),
  Top_NTL_CellType   = nt_comp_cols[max.col(mean_cpm_by_comp[all_sig_genes, nt_comp_cols, drop = FALSE], ties.method = "first")],
  Top_NTL_Cell_CPM   = apply(mean_cpm_by_comp[all_sig_genes, nt_comp_cols, drop = FALSE], 1, max),
  Top_NTL_Cell_Pct   = sapply(seq_along(all_sig_genes), function(i) {
    pct_by_comp[all_sig_genes[i], nt_comp_cols[which.max(mean_cpm_by_comp[all_sig_genes[i], nt_comp_cols])]]
  })
) %>%
  mutate(
    # 定义“污染易感基因”: 在癌旁整体表达显著高于肿瘤 (Log2FC > 0.8 且癌旁 CPM > 10)
    Admixture_Prone = (Log2FC_NTL_vs_T > 0.8 & PB_FullNTL_CPM > 10)
  )

# 把每个基因映射回所属免疫特征
sig_gene_map <- map_dfr(names(sig_list), function(s) {
  tibble(Signature = s, Gene = intersect(sig_list[[s]], all_sig_genes))
}) %>%
  left_join(gene_source_df, by = "Gene")

# 打印大家最关心的两个高敏感特征里的 Top 污染易感基因
for (target_sig in c("T_cell_exhaustion", "HM_Complement", "HM_TNFA_via_NFKB")) {
  cat(sprintf("\n--- [%s] 排名前 10 的癌旁高本底“污染易感基因” ---\n", target_sig))
  top_g <- sig_gene_map %>%
    filter(Signature == target_sig) %>%
    arrange(desc(Log2FC_NTL_vs_T)) %>%
    head(10) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    select(Gene, PB_Tumor_CPM, PB_FullNTL_CPM, Log2FC_NTL_vs_T, Top_NTL_CellType, Top_NTL_Cell_Pct)
  print(as.data.frame(top_g), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 5. 构建零循环论证的“正常肝实质混入评分 (NHAS, 30 基因)”并在三大 Bulk 队列验收
# ------------------------------------------------------------------------------
message("\n>>> [5/6] 正在从单细胞构建无免疫循环论证的 NHAS (Normal Hepatocyte Admixture Score)...")

idx_N_hep   <- which(meta_tn$Site %in% unique(meta$Site[is_N]) & grepl("Hep", meta_tn$CellType, ignore.case = TRUE))
idx_T_hep   <- which(meta_tn$Site %in% unique(meta$Site[is_T]) & grepl("Hep", meta_tn$CellType, ignore.case = TRUE))
idx_Non_hep <- which(!grepl("Hep", meta_tn$CellType, ignore.case = TRUE))

m_N_hep   <- Matrix::rowMeans(sp_tn_cpm[, idx_N_hep,   drop = FALSE])
m_T_hep   <- Matrix::rowMeans(sp_tn_cpm[, idx_T_hep,   drop = FALSE])
m_Non_hep <- Matrix::rowMeans(sp_tn_cpm[, idx_Non_hep, drop = FALSE])

p_N_hep   <- Matrix::rowMeans(sp_tn[, idx_N_hep,   drop = FALSE] > 0)
p_Non_hep <- Matrix::rowMeans(sp_tn[, idx_Non_hep, drop = FALSE] > 0)

# 确保候选基因在三大 Bulk 队列中均存在，且与 24 个免疫特征严格零交集
tcga_obj  <- readRDS("01_clean_data/tcga_lihc_clean.rds")
h18_obj   <- readRDS("01_clean_data/hccdb18_clean.rds")
gse_obj   <- readRDS("01_clean_data/gse14520_clean.rds")
common_3bulk <- Reduce(intersect, list(rownames(tcga_obj$T_tpm), rownames(h18_obj$T_tpm), rownames(gse_obj$T_lin)))
immune_union <- unique(unlist(sig_list))

nhas_cand <- tibble(
  Gene                 = rownames(sp_tn),
  Mean_NormHep_CPM     = m_N_hep,
  Mean_TumorHep_CPM    = m_T_hep,
  Mean_NonHep_CPM      = m_Non_hep,
  Pct_NormHep          = p_N_hep,
  Pct_NonHep           = p_Non_hep,
  Log2FC_vs_NonHep     = log2((m_N_hep + 1) / (m_Non_hep + 1)),
  Log2FC_vs_TumorHep   = log2((m_N_hep + 1) / (m_T_hep + 1))
) %>%
  filter(
    Gene %in% common_3bulk,
    !(Gene %in% immune_union),      # 强制与所有免疫特征零交集！
    !grepl("^MT-|^RP[SL]", Gene),   # 剔除线粒体与核糖体基因
    Pct_NormHep >= 0.40,            # 在正常肝细胞中稳定检出
    Pct_NonHep  <= 0.05,            # 在所有免疫/内皮/成纤维细胞中近乎零表达
    Log2FC_vs_NonHep   >= 3.0,      # 相比非肝实质高度特异
    Log2FC_vs_TumorHep >= 1.0       # 在正常肝细胞中显著高于恶性肝癌细胞
  ) %>%
  mutate(Specificity_Score = Log2FC_vs_TumorHep * Pct_NormHep) %>%
  arrange(desc(Specificity_Score))

nhas_genes <- head(nhas_cand$Gene, 30)
stopifnot("NHAS 基因与免疫特征存在交集！" = length(intersect(nhas_genes, immune_union)) == 0)

message(sprintf("    成功筛选出 %d 个零循环论证 NHAS 基因 (与 24 个免疫特征交集 = 0):", length(nhas_genes)))
cat("   ", paste(nhas_genes, collapse = ", "), "\n")

# 在三大 Bulk 模拟队列中检验 NHAS 对真实混入比例 alpha (0%–50%) 的单调追踪精度
validate_nhas <- function(T_mat, N_mat, label) {
  a_grid <- c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
  # 对每位患者拼接 7 个梯度样本计算单样本 ssGSEA，检验患者内与队列整体相关性
  nhas_mat <- sapply(a_grid, function(a) {
    m <- log2((1 - a) * T_mat + a * N_mat + 1)
    p <- ssgseaParam(as.matrix(m), list(NHAS = nhas_genes), minSize = 5, normalize = FALSE)
    as.numeric(gsva(p, verbose = FALSE)[1, ])
  })
  # 每位患者内部 NHAS 随 alpha 递增的 Spearman 相关系数中位数
  within_rhos <- apply(nhas_mat, 1, function(v) cor(a_grid, v, method = "spearman"))
  # 队列层面 NHAS 差值 (alpha - 0) 与 alpha 的总体相关性
  delta_vec <- as.numeric(sweep(nhas_mat, 1, nhas_mat[, 1], "-"))
  alpha_vec <- rep(a_grid * 100, each = nrow(nhas_mat))
  overall_r <- cor(alpha_vec, delta_vec, method = "spearman")
  tibble(
    Cohort = label,
    Median_WithinPatient_Rho = round(median(within_rhos, na.rm = TRUE), 3),
    Pct_Patients_Rho_gt_0.90 = round(mean(within_rhos >= 0.90, na.rm = TRUE) * 100, 1),
    Overall_DeltaNHAS_Rho    = round(overall_r, 3)
  )
}

nhas_val_tbl <- bind_rows(
  validate_nhas(tcga_obj$T_tpm, tcga_obj$N_tpm, "TCGA-LIHC (n=50)"),
  validate_nhas(h18_obj$T_tpm,  h18_obj$N_tpm,  "ICGC-LIRI-JP (n=175)"),
  validate_nhas(gse_obj$T_lin,  gse_obj$N_lin,  "GSE14520 (n=214)")
)
message("\n    [NHAS 量尺金标准验收结果]:")
print(as.data.frame(nhas_val_tbl), row.names = FALSE)

# 保存第 3 关全部核心成果
saveRDS(
  list(
    paired_pts       = paired_pts,
    pb_res           = pb_res,
    decomp_30        = decomp_30,
    sig_gene_map     = sig_gene_map,
    mean_cpm_by_comp = mean_cpm_by_comp,
    pct_by_comp      = pct_by_comp,
    nhas_genes       = nhas_genes,
    nhas_cand        = nhas_cand,
    nhas_val_tbl     = nhas_val_tbl
  ),
  file = "01_clean_data/03_scrna_mechanism_and_nhas.rds"
)
write.csv(decomp_30,    "03_results/03_scrna_mechanism_decomposition.csv", row.names = FALSE)
write.csv(sig_gene_map, "03_results/03_signature_gene_source_attribution.csv", row.names = FALSE)
write.csv(nhas_cand,    "03_results/03_nhas_candidate_genes.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. 绘制 Figure 3: 单细胞机制拆解图 (Arm A vs Arm B 哑铃图 + 成分堆叠条形图)
# ------------------------------------------------------------------------------
message("\n>>> [6/6] 正在绘制 Figure 3...")

sig_ord3 <- rev(decomp_30$Signature)

p3a <- pb_res %>%
  filter(Method == "ssGSEA", Alpha_Pct == 30) %>%
  mutate(Signature = factor(Signature, levels = sig_ord3)) %>%
  ggplot(aes(x = Std_Drift, y = Signature)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(aes(group = Signature), color = "grey60", linewidth = 1) +
  geom_point(aes(color = Arm), size = 3.2) +
  scale_color_manual(
    values = c("ArmA_Pure_Hepatocytes" = "#1F77B4", "ArmB_Full_Adjacent_Liver" = "#D62728"),
    labels = c("Arm A: Pure Normal Hepatocytes (30%)", "Arm B: Full Adjacent Liver (30%)"),
    name   = "Pseudo-bulk Admixture Arm"
  ) +
  labs(
    title    = "A. Pure Hepatocytes vs. Full Adjacent Liver",
    subtitle = sprintf("Paired scRNA-seq pseudo-bulk (GSE149614, n=%d patients)", length(paired_pts)),
    x        = "Standardized Score Drift (ssGSEA, 30%)",
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# B 图: 展示每个高敏感特征中“污染易感基因”的癌旁主导细胞类型占比
prone_summary <- sig_gene_map %>%
  filter(Admixture_Prone) %>%
  mutate(Cell_Group = sub("^N_", "", Top_NTL_CellType)) %>%
  count(Signature, Cell_Group) %>%
  mutate(Signature = factor(Signature, levels = sig_ord3))

p3b <- ggplot(prone_summary, aes(x = n, y = Signature, fill = Cell_Group)) +
  geom_col(width = 0.72) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_brewer(palette = "Set2", name = "Dominant NTL Cell Source") +
  labs(
    title    = "B. Cellular Origin of Admixture-Prone Genes",
    subtitle = "Signature genes enriched in adjacent liver (log2FC > 0.8, CPM > 10)",
    x        = "Number of Admixture-Prone Genes",
    y        = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

p3 <- (p3a | p3b) + plot_layout(widths = c(1.3, 1))

ggsave("04_figures/Fig3_scRNA_Mechanism_Decomposition.pdf", p3, width = 13.5, height = 8)
ggsave("04_figures/Fig3_scRNA_Mechanism_Decomposition.png", p3, width = 13.5, height = 8, dpi = 300)

message(">>> 全部完成！请查看 04_figures/Fig3_scRNA_Mechanism_Decomposition.png。")
