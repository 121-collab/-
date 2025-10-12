# ============================================================
# Liao et al., 2020 人肾脏 scRNA-seq 复现 — 单文件脚本
# 关键参数与阈值对齐论文：QC、Harmony(20 PCs)、res=0.25/0.8
# 参考：Fig.1–4, Table 2:contentReference[oaicite:1]{index=1}
# ============================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix)
  library(harmony)              # 若未装：remotes::install_github("immunogenomics/harmony")
  library(ggplot2); library(patchwork); library(pheatmap); library(scales)
})

# ------------------ 可调参数 ------------------
set.seed(123)
DATA_DIRS <- c(
  "D:/scRNAseq/GSE131685_RAW/kidney1",
  "D:/scRNAseq/GSE131685_RAW/kidney2",
  "D:/scRNAseq/GSE131685_RAW/kidney3"
)
N_HVG            <- 2000        # 可变基因数量（1500/1000 也可，越小越省内存）
HARMONY_DIMS     <- 1:20        # 论文使用 20 PCs:contentReference[oaicite:2]{index=2}
RES_MAIN         <- 0.25        # 主聚类分辨率（论文 0.25）
RES_FINE         <- 0.8         # 细分分辨率（用于 NKT/T 等）
TARGET_DS        <- 10000       # 可视化降采样目标细胞数（不影响计算）
OUT_DIR          <- "kidney_repro_out"
dir.create(OUT_DIR, showWarnings = FALSE)
# ---------------------------------------------

# =============== 1) 读入与合并（立刻 JoinLayers） ===============
mats <- lapply(DATA_DIRS, Read10X)
objs <- Map(function(m, i) CreateSeuratObject(m, project = paste0("Kidney", i)), mats, seq_along(mats))
kidney_merged <- Reduce(function(a,b) merge(a, b, add.cell.ids=c(a@project.name, b@project.name)), objs)

DefaultAssay(kidney_merged) <- "RNA"
# 关键：Seurat v5 立刻合并 layer，避免 counts.Kidney1/… 残留
kidney_merged <- JoinLayers(kidney_merged, assay="RNA")
stopifnot(!any(grepl("\\.", Layers(kidney_merged[["RNA"]]))))  # 现在应只剩 "counts"

# =============== 2) 论文阈值 QC（Fig.2a–b） ===============
mt_pattern <- if (any(grepl("^MT-", rownames(kidney_merged)))) "^MT-" else "^mt-"
kidney_merged[["percent.mt"]] <- PercentageFeatureSet(kidney_merged, pattern = mt_pattern)

k <- subset(
  kidney_merged,
  subset = nFeature_RNA >= 200 & nFeature_RNA <= 2500 & percent.mt <= 30   # 与论文一致:contentReference[oaicite:3]{index=3}
)
k$batch <- k$orig.ident

# QC 可视化（与 Fig.2a–b 对照）
p_qc1 <- FeatureScatter(k, feature1="nCount_RNA", feature2="percent.mt", group.by="batch") +
  geom_hline(yintercept = 30, linetype="dashed") + ggtitle("percent.mt vs nCount (cut=30%)")
p_qc2 <- FeatureScatter(k, feature1="nCount_RNA", feature2="nFeature_RNA", group.by="batch") +
  ggtitle("nFeature vs nCount")
p_vln <- VlnPlot(k, features=c("nFeature_RNA","nCount_RNA","percent.mt"), group.by="batch", pt.size=0, ncol=3)
print((p_qc1 | p_qc2) / p_vln)

# =============== 3) Normalize → HVG → Scale(HVG) → PCA ===============
k <- NormalizeData(k)
k <- FindVariableFeatures(k, selection.method="vst", nfeatures=N_HVG)
k <- ScaleData(k, features = VariableFeatures(k), verbose = FALSE)
k <- RunPCA(k, features = VariableFeatures(k), npcs = 50, verbose = FALSE)

# =============== 4) Harmony 批次校正（20 PCs）:contentReference[oaicite:4]{index=4} ===============
rh_formals <- names(formals(getS3method("RunHarmony","Seurat")))
if ("dims" %in% rh_formals) {
  k <- RunHarmony(
    object = k, group.by.vars = "batch",
    reduction = "pca", dims = HARMONY_DIMS,
    assay.use = DefaultAssay(k),
    plot_convergence = FALSE, verbose = TRUE
  )
} else {
  k <- RunHarmony(
    object = k, group.by.vars = "batch",
    reduction.use = "pca", dims.use = HARMONY_DIMS,
    assay.use = DefaultAssay(k),
    plot_convergence = FALSE, verbose = TRUE
  )
}
print(DimPlot(k, reduction="harmony", group.by="batch") + ggtitle("Batch mixing in Harmony space"))

# =============== 5) UMAP/邻域/聚类（基于 Harmony 维度） ===============
k <- RunUMAP(k, reduction="harmony", dims=HARMONY_DIMS)
k <- FindNeighbors(k, reduction="harmony", dims=HARMONY_DIMS)
k <- FindClusters(k, resolution=RES_MAIN)

umap_main <- DimPlot(k, reduction="umap", label=TRUE, repel=TRUE) +
  ggtitle(sprintf("UMAP (res=%.2f) — expect ~10 clusters", RES_MAIN))
print(umap_main)

# =============== 6) 簇×批次分层降采样（仅用于可视化，不影响统计） ===============
ids_keep <- unlist(lapply(levels(Idents(k)), function(cl){
  cells_cl <- WhichCells(k, idents = cl)
  if (length(cells_cl) == 0) return(character(0))
  n_cl <- round(TARGET_DS * length(cells_cl) / ncol(k))
  if (n_cl >= length(cells_cl)) return(cells_cl)
  byb <- split(cells_cl, k$batch[cells_cl])
  prop <- sapply(byb, length) / length(cells_cl); take <- pmax(1, round(n_cl * prop))
  unlist(mapply(function(v,n) sample(v, size=min(n, length(v))), byb, take, SIMPLIFY=FALSE), use.names=FALSE)
}))
k_ds <- subset(k, cells = unique(ids_keep))
print(DimPlot(k_ds, reduction="umap", label=TRUE) + ggtitle(sprintf("UMAP downsampled (n=%d)", ncol(k_ds))))

# =============== 7) 文献 marker 集 & 自动注释（Fig.1d 对齐） ===============
DefaultAssay(k_ds) <- "RNA"
markers_ref <- list(
  "Proximal convoluted tubule"        = c("SLC34A1","LRP2","ALDOB","SLC5A2","CUBN","GPX3"),
  "Proximal tubule"                   = c("SLC34A1","ALDOB","LRP2","DCXR","GPX3"),
  "Proximal straight tubule"          = c("SLC13A3","SLC22A8","SLC22A7"),
  "NK-T cells"                        = c("GNLY","NKG7","PRF1","GZMB","CD3D","CD3E","IL7R"),
  "Monocytes"                         = c("LST1","S100A8","S100A9","FCN1","LYZ","MS4A7"),
  "Glomerular parietal epithelial cells" = c("KRT8","KRT18","KRT19","EPCAM"),
  "Distal tubule cells"               = c("SLC12A3","TRPM6","CLDN16"),
  "Collecting duct principal cells"   = c("AQP2","SCNN1G","FXYD4"),
  "Collecting duct intercalated cells"= c("ATP6V1B1","ATP6V0D2","ATP6V1G3","FOXI1"),
  "B cells"                           = c("MS4A1","CD79A","CD79B","CD74")
)
markers_ref <- lapply(markers_ref, \(v) intersect(v, rownames(k_ds)))

avg <- AverageExpression(k_ds, assays="RNA", slot="data", group.by="seurat_clusters")$RNA
colnames(avg) <- sub("^g", "", colnames(avg))                 # v5 可能有 g 前缀
avg_z <- t(scale(t(avg)))                                     # 行 z-score（基因内）

score_mat <- sapply(names(markers_ref), function(tp){
  genes <- markers_ref[[tp]]; if (length(genes)==0) return(rep(NA_real_, ncol(avg_z)))
  colMeans(avg_z[genes, , drop = FALSE], na.rm = TRUE)
})
rownames(score_mat) <- colnames(avg_z)
best_type <- colnames(score_mat)[apply(score_mat, 1, which.max)]    # 每簇最佳类型
names(best_type) <- rownames(score_mat)

# 写回注释
# ---------- 安全写入 annot_ref（替换原来的几行） ----------
# 前置：已得到每个簇的最佳类型映射 best_type，形如 c("0"="Proximal tubule", "1"="PCT", ...)

# 1) 先把当前簇号写到 meta 里（按细胞；字符型）
k_ds$cluster_tmp <- as.character(Idents(k_ds))   # 每个细胞对应的簇号，如 "0","1",...

# 2) 用 best_type 做“逐细胞映射”得到 annot_ref（长度 = 细胞数）
annot_vec <- unname(best_type[k_ds$cluster_tmp])  # 向量长度与细胞数一致
annot_vec[is.na(annot_vec)] <- "Unknown"          # 没匹配到的填 Unknown

# 3) VERY IMPORTANT: 给向量设置细胞名，确保与对象严格对齐
names(annot_vec) <- colnames(k_ds)

# 4) 防御式检查：重叠数量、前几项预览
overlap_n <- sum(names(annot_vec) %in% colnames(k_ds))
if (overlap_n == 0) {
  stop("annot_vec 与对象的细胞名完全不重叠：请检查是否对同一个对象进行了子集/重聚类但没同步生成映射。")
}
message("将写入 annot_ref：长度=", length(annot_vec), "；细胞名重叠=", overlap_n)

# 5) 两种安全写法，任选其一（优先 A）
# A) 向量法（要求 names(annot_vec) 为细胞名）
k_ds <- AddMetaData(k_ds, metadata = annot_vec, col.name = "annot_ref")

# B) data.frame 法（与 rownames 对齐；若你更放心用表格，可启用此法并注释掉 A）
# df_meta <- data.frame(annot_ref = annot_vec, row.names = names(annot_vec), check.names = FALSE)
# k_ds <- AddMetaData(k_ds, metadata = df_meta)

# 6) 验证写入成功
print(head(k_ds@meta.data[, c("cluster_tmp", "annot_ref")], 3))
print(table(k_ds$annot_ref))


print(table(k_ds$annot_ref))

# UMAP（按文献注释上色）
pal <- c(
  "Proximal convoluted tubule"="#E69F00","Proximal tubule"="#F28E2B",
  "Proximal straight tubule"="#2CA02C","NK-T cells"="#E15759","Monocytes"="#FF9DA7",
  "Glomerular parietal epithelial cells"="#8DD3C7","Distal tubule cells"="#66C2A5",
  "Collecting duct principal cells"="#A6CEE3","Collecting duct intercalated cells"="#6A3D9A",
  "B cells"="#1F78B4","Unknown"="grey70"
)
p_umap_annot <- DimPlot(k_ds, reduction="umap", group.by="annot_ref", label=TRUE, repel=TRUE) +
  scale_color_manual(values=pal, na.value="grey70") +
  labs(title="UMAP — literature-based annotation", color=NULL)
print(p_umap_annot)

# 饼图（类型占比）
df_pie <- as.data.frame(table(k_ds$annot_ref)); colnames(df_pie) <- c("type","n")
df_pie$frac <- df_pie$n / sum(df_pie$n)
p_pie <- ggplot(df_pie, aes(x="", y=frac, fill=type)) +
  geom_col() + coord_polar("y") + scale_fill_manual(values=pal) +
  theme_void() + labs(title="Kidney cell-type proportions")
print(p_pie)

# 文献风格 marker 热图（群平均 + 基因内 z；颜色截断）
markers_trim <- lapply(markers_ref, \(v) head(v, 8))
features_ordered <- unique(unlist(markers_trim)); features_ordered <- intersect(features_ordered, rownames(k_ds))
k_ds <- ScaleData(k_ds, features=features_ordered, verbose=FALSE)
ht_ref <- DoHeatmap(k_ds, features=features_ordered, group.by="seurat_clusters",
                    raster=TRUE, slot="scale.data", draw.lines=FALSE) +
  scale_fill_gradientn(colors=c("navy","black","yellow")) +
  ggtitle("Literature markers heatmap (cluster means, z-scored)")
print(ht_ref)



# =============== 8) 集合管主/闰验证（Fig.4a） ===============
FeaturePlot(k_ds, features = c("AQP2","ATP6V1B1","ATP6V0D2","ATP6V1G3"), ncol=2)

# =============== 9) NKT vs T 细分（Fig.4b–g） ===============
k_ds <- FindClusters(k_ds, resolution = RES_FINE)
genes_nkt <- intersect(c("GNLY","NKG7","GZMB","PRF1"), rownames(k_ds))
genes_t   <- intersect(c("CD3D","CD3E","IL7R","CCR7"), rownames(k_ds))

k_ds <- AddModuleScore(k_ds, features=list(genes_nkt), name="NKT_Score", assay="RNA")
k_ds <- AddModuleScore(k_ds, features=list(genes_t),   name="T_Score",   assay="RNA")

delta <- 0.1
nkt   <- k_ds$NKT_Score1; tlin <- k_ds$T_Score1
k_ds$NT_class <- factor(ifelse(nkt > tlin + delta, "NKT cells",
                               ifelse(tlin > nkt + delta, "T cells", "Other")),
                        levels=c("NKT cells","T cells","Other"))
nt <- subset(k_ds, subset = NT_class %in% c("NKT cells","T cells"))
bal <- min(table(nt$NT_class)); set.seed(123)
keep <- unlist(tapply(colnames(nt), nt$NT_class, function(v) sample(v, bal)))
nt_bal <- subset(nt, cells = keep)

p_umap_nt <- DimPlot(nt_bal, reduction="umap", group.by="NT_class") +
  scale_color_manual(values=c("NKT cells"="#E76F51", "T cells"="#2A9D8F")) +
  labs(title="NKT vs T on UMAP", color=NULL)
print(p_umap_nt)

feat_violin <- intersect(c("CD3D","CD3E","GNLY","NKG7","IL7R"), rownames(nt_bal))
ymax <- max(FetchData(nt_bal, vars = feat_violin))
vp <- VlnPlot(nt_bal, features=feat_violin, group.by="NT_class",
              pt.size=0.6, combine=FALSE, slot="data")
vp <- lapply(vp, function(p) p + coord_cartesian(ylim=c(0, ymax)) +
               xlab("") + ylab("Expression Level") + theme(legend.position="none"))
print(wrap_plots(vp, ncol=3))

# =============== 10) PT 三亚群（Fig.3a 的前置验证） ===============
pt_cells <- WhichCells(k_ds, expression = SLC34A1 > 0 | LRP2 > 0 | ALDOB > 0)
pt <- subset(k_ds, cells = pt_cells)
pt <- FindNeighbors(pt, reduction = "harmony", dims = HARMONY_DIMS)
pt <- FindClusters(pt, resolution = 0.3)
print(DimPlot(pt, reduction = "umap", label=TRUE) + ggtitle("PT subclusters"))
FeaturePlot(pt, features=c("SLC5A2","LRP2","SLC22A8","ALDOB","GPX3"), ncol=3)

# =============== 11) 论文编号（1–10）上图 & 饼图（对齐 Fig.1） ===============
paper_order <- c(
  "Proximal convoluted tubule"         = 1,
  "Proximal tubule"                    = 2,
  "Proximal straight tubule"           = 3,
  "NK-T cells"                         = 4,
  "Monocytes"                          = 5,
  "Glomerular parietal epithelial cells" = 6,
  "Distal tubule cells"                = 7,
  "Collecting duct principal cells"    = 8,
  "B cells"                            = 9,
  "Collecting duct intercalated cells" = 10
)
num_vec <- setNames(paper_order[as.character(k_ds$annot_ref)], colnames(k_ds))
k_ds$cluster_paper <- factor(num_vec, levels = as.character(1:10))
pal_by_type <- c(
  "Proximal convoluted tubule"="#E69F00",
  "Proximal tubule"="#F28E2B",
  "Proximal straight tubule"="#2CA02C",
  "NK-T cells"="#E15759",
  "Monocytes"="#FF9DA7",
  "Glomerular parietal epithelial cells"="#8DD3C7",
  "Distal tubule cells"="#66C2A5",
  "Collecting duct principal cells"="#A6CEE3",
  "B cells"="#1F78B4",
  "Collecting duct intercalated cells"="#6A3D9A"
)
pal_paper <- setNames(pal_by_type[names(paper_order)], as.character(paper_order))

Idents(k_ds) <- "cluster_paper"
p_umap_paper <- DimPlot(k_ds, reduction="umap", label=TRUE, label.size=5, repel=TRUE) +
  scale_color_manual(values = pal_paper, drop = FALSE) +
  labs(title = "UMAP — paper-style clusters (1–10)", color = NULL) +
  theme_bw(base_size = 12) + theme(panel.grid = element_blank())
print(p_umap_paper)

pie_df <- as.data.frame(table(k_ds$cluster_paper), stringsAsFactors = FALSE)
names(pie_df) <- c("cluster_id","n")
pie_df <- pie_df %>% filter(!is.na(cluster_id)) %>%
  mutate(frac = n / sum(n))
p_pie_paper <- ggplot(pie_df, aes(x="", y=frac, fill=cluster_id)) +
  geom_col(width=1) + coord_polar("y") +
  scale_fill_manual(values = pal_paper, drop = FALSE) +
  theme_void(base_size = 11) + labs(title = "Kidney cell types (paper-style 1–10)")
print(p_pie_paper)
## 📚 参考资料
- Liao et al., *Nature Medicine*, 2020  
- GEO: [GSE131685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131685)  
- 10x Genomics 官方文档  
- [Seurat 官方教程](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
