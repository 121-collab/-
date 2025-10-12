# 🧬 Liao et al., 2020 — Human Kidney scRNA-seq Reproduction

**论文信息**  
> Liao et al., *Nature Medicine*, 2020  
> DOI: [10.1038/s41591-020-0818-9](https://doi.org/10.1038/s41591-020-0818-9)  
> 数据集：GEO [GSE131685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131685)

---

## 1️⃣ 数据信息

| 项目 | 内容 |
|------|------|
| 平台 | 10x Genomics Chromium v2 |
| 物种 | *Homo sapiens* |
| 数据类型 | scRNA-seq |
| 原始文件 | FASTQ (SRA: SRR10377488–90) |
| 参考基因组 | GRCh38 |
| 分析工具 | Cell Ranger, Seurat |

---

## 2️⃣ 数据下载与预处理

```bash
# 创建目录
mkdir -p ~/projects/Liao2020_kidney
cd ~/projects/Liao2020_kidney

# 写入 SRA 号
cat > SraAccList.txt <<EOF
SRR10377488
SRR10377489
SRR10377490
EOF

# 下载与转换
prefetch --option-file SraAccList.txt
fasterq-dump --split-files SRR10377488 SRR10377489 SRR10377490
```

---

## 3️⃣ Cell Ranger 分析

```bash
# 构建参考基因组
cellranger mkref --genome=GRCh38 --fasta=GRCh38.fa --genes=genes.gtf

# 运行 count
cellranger count   --id=Liao_kidney   --transcriptome=/path/to/refdata-GRCh38-2020A   --fastqs=/path/to/fastq   --sample=SRR10377488,SRR10377489,SRR10377490   --localcores=16 --localmem=128
```

输出：
- `outs/filtered_feature_bc_matrix/`
- `outs/raw_feature_bc_matrix/`
- `outs/web_summary.html`
- `outs/metrics_summary.csv`

---

## 4️⃣ R 分析 (Seurat)

```R
library(Seurat)
library(dplyr)

# 读取数据
kidney <- Read10X(data.dir = "filtered_feature_bc_matrix/")
kidney <- CreateSeuratObject(counts = kidney, project = "Liao2020_kidney")

# 质控
kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^MT-")
kidney <- subset(kidney, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# 标准化与高变基因
kidney <- NormalizeData(kidney)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)

# PCA + 聚类 + 可视化
kidney <- ScaleData(kidney)
kidney <- RunPCA(kidney)
kidney <- FindNeighbors(kidney, dims = 1:20)
kidney <- FindClusters(kidney, resolution = 0.5)
kidney <- RunUMAP(kidney, dims = 1:20)
DimPlot(kidney, reduction = "umap", label = TRUE)

# 差异基因分析
markers <- FindAllMarkers(kidney, only.pos = TRUE)
write.csv(markers, "all_markers.csv", row.names = FALSE)
```

---

## 5️⃣ 输出验证
- 对照论文中的 QC 图与 cluster 数量（Fig.1–3）
- 验证主要细胞类型（PT, TAL, DCT, CD, Podocyte 等）
- 生成 `metrics_summary.csv` 并记录核心指标

---

## 📚 参考资料
- Liao et al., *Nature Medicine*, 2020  
- GEO: [GSE131685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131685)  
- 10x Genomics 官方文档  
- [Seurat 官方教程](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
