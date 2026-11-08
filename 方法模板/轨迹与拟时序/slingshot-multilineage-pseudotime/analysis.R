## Slingshot 单细胞指定起始群多谱系拟时序推断
library(Seurat)
library(SingleCellExperiment)
library(slingshot)

sce_seurat <- readRDS("Data/bcell_seurat.rds")
Idents(sce_seurat) <- "B_type"

seuratObject <- as.SingleCellExperiment(sce_seurat)
seuratObject <- slingshot(
  seuratObject,
  clusterLabels = "B_type",
  reducedDim = "TSNE",
  start.clus = "Naive_B",
  stretch = 0,
  approx_points = 100
)

psts <- as.data.frame(slingPseudotime(seuratObject))

# 每细胞在其所属谱系上的伪时间取最大值（未分配谱系记 0），得到跨谱系可比的 pseudotime0
seuratObject$pseudotime0 <- apply(psts, 1, function(x) max(x, 0, na.rm = TRUE))

saveRDS(seuratObject, "results/slingshot_sce.rds")
write.csv(psts, "results/slingshot_pseudotime.csv")
