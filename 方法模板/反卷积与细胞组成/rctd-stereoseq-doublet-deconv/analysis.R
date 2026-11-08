## spacexr (RCTD) 高分辨率空间组（Stereo-seq/Slide-seq）doublet 模式反卷积
library(spacexr)
library(Seurat)
library(dplyr)

snRNA <- readRDS("Data/snRNA_seurat.rds")
ST <- readRDS("Data/stereo_seurat.rds")

# 单细胞核参考提供细胞类型表达谱先验
reference <- Reference(
  GetAssayData(snRNA, assay = "RNA", layer = "counts"),
  cell_types = snRNA$subclass,
  nUMI = snRNA$nCount_RNA
)

coords <- ST@meta.data[, c("x", "y")]
countsSP <- GetAssayData(ST, layer = "counts")
umiSP <- Matrix::colSums(countsSP)
puck <- SpatialRNA(coords, countsSP, nUMI = umiSP)

# 将参考表达谱限制在与空间数据共有的基因上，避免特征空间错位
reference@counts <- reference@counts[rownames(reference@counts) %in% rownames(countsSP), ]
reference@nUMI <- Matrix::colSums(reference@counts)

myRCTD <- create.RCTD(puck, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

# doublet 模式为每个 bin 给出 1-2 个细胞类型（first/second）与 singlet/doublet 分类
saveRDS(myRCTD@results, "results/rctd_doublet_results.rds")
