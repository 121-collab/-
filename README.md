# 🔬 Bioinformatics Reproduction Notes

> 个人科研复现笔记库 —— 专注于 Illumina 二代测序与单细胞转录组（scRNA-seq）数据的复现与分析  
> Author: 银河 | Role: 生信初学者 / scRNA-seq 技术学习者  
> 环境：Linux 服务器 + Cell Ranger + R + Seurat

---

## 📁 仓库结构

```
bio-replication-notes/
│
├── README.md                    # 总览文档（本文件）
│
├── Liao2020_kidney/             # 当前复现项目：Liao et al., Nature Medicine 2020
│   ├── README.md                # Liao2020 项目复现笔记（完整流程）
│   ├── run_cellranger.sh        # 数据分析脚本模板
│   ├── analysis_seurat.R        # Seurat 分析主脚本
│   └── results/                 # 输出结果（质控、聚类、图表等）
│
└── future_projects/             # 预留后续复现项目
    ├── GSE123516/
    ├── PBMC10x/
    └── ...
```

---

## 🧩 当前复现项目

### 📖 Liao et al., 2020 — *Single-cell landscape of human kidney under healthy and diseased conditions*
- **期刊**: Nature Medicine (2020)  
- **DOI**: [10.1038/s41591-020-0818-9](https://doi.org/10.1038/s41591-020-0818-9)  
- **数据来源**: GEO [GSE131685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131685)  
- **SRA 号**: SRR10377488, SRR10377489, SRR10377490  
- **物种**: *Homo sapiens*  
- **平台**: 10x Genomics Chromium v2  
- **分析内容**: Cell Ranger + Seurat 全流程复现  

📘 详细复现步骤请见 → [Liao2020_kidney/README.md](./Liao2020_kidney/README.md)

---

## 🚀 使用说明

### ▶️ 快速启动

```bash
# 克隆仓库
git clone https://github.com/<your-username>/bio-replication-notes.git
cd bio-replication-notes/Liao2020_kidney

# 运行 Cell Ranger 分析
bash run_cellranger.sh

# 运行 R 分析
Rscript analysis_seurat.R
```

---

## 🧠 未来计划

| 项目 | 期刊 | 数据集 | 状态 |
|------|------|--------|------|
| Liao2020_kidney | Nature Medicine 2020 | GSE131685 | ✅ 已完成 |
| GSE123516_intestine | Nature 2019 | GSE123516 | 🔄 计划中 |
| Zheng2017_PBMC | Cell 2017 | GSE99254 | 🔄 计划中 |

---

## 💡 学习目标

- 从零开始掌握 **scRNA-seq 全流程分析**
- 熟练使用 **Linux + Cell Ranger + Seurat**
- 独立编写与理解每个分析脚本
- 实现可复现、结构化的科研记录

---

## 🧾 License

本仓库仅用于科研学习与方法复现，不含原始数据。  
请遵守对应论文与数据库的公开数据使用协议。
