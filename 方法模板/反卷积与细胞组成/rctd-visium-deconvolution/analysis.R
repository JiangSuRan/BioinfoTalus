## spacexr (RCTD) 10x Visium 空间转录组细胞类型反卷积
library(spacexr)
library(Seurat)

set.seed(42)

seu_ref <- readRDS("Data/sc_reference_seurat.rds")
seu_sp <- readRDS("Data/visium_seurat.rds")

# 单细胞参考提供各细胞类型的表达谱先验
reference <- Reference(
  GetAssayData(seu_ref, assay = "RNA", layer = "counts"),
  cell_types = seu_ref$cell_type,
  nUMI = seu_ref$nCount_RNA
)

# SpatialRNA 由 spot 坐标与计数构成，坐标顺序须与计数列一致
coords <- GetTissueCoordinates(seu_sp, cols = c("imagecol", "imagerow"), scale = NULL)
spatialRNA <- SpatialRNA(
  coords,
  GetAssayData(seu_sp, assay = "Spatial", layer = "counts"),
  nUMI = seu_sp$nCount_Spatial
)

myRCTD <- create.RCTD(spatialRNA, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

saveRDS(myRCTD, "results/rctd_visium.rds")

# 权重按 spot 行归一化后导出，行=spot、列=细胞类型，行为该 spot 的细胞类型组成
weights <- as.data.frame(as.matrix(myRCTD@results$weights))
norm_weights <- sweep(weights, 1, rowSums(weights), "/")
write.csv(norm_weights, "results/rctd_weights.csv")
