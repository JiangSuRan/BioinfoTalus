## Visium HD 高分辨率多尺度分箱加载、LeverageScore Sketch 降采样映射与 Banksy 空间域聚类
library(Seurat)
library(SeuratWrappers)
library(Banksy)

# 同时载入 8 um 亚细胞级与 16 um 增强覆盖率分箱矩阵，以 LeverageScore Sketch 加速全图谱聚类并投影回全量斑点
hd_obj <- Load10X_Spatial(data.dir = "Data/visium_hd_out", bin.size = c(8, 16))
DefaultAssay(hd_obj) <- "Spatial.008um"

hd_obj <- NormalizeData(hd_obj)
hd_obj <- FindVariableFeatures(hd_obj)
hd_obj <- ScaleData(hd_obj)

hd_obj <- SketchData(object = hd_obj, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(hd_obj) <- "sketch"
hd_obj <- FindVariableFeatures(hd_obj)
hd_obj <- ScaleData(hd_obj)
hd_obj <- RunPCA(hd_obj, assay = "sketch", reduction.name = "pca.sketch")
hd_obj <- FindNeighbors(hd_obj, assay = "sketch", reduction = "pca.sketch", dims = 1:30)
hd_obj <- FindClusters(hd_obj, cluster.name = "seurat_cluster.sketched", resolution = 1.5)
hd_obj <- RunUMAP(hd_obj, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:30)

hd_obj <- ProjectData(
  object = hd_obj,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# 联合斑点自身表达与邻域微环境表达（lambda = 0.2）区分精细空间组织结构域
DefaultAssay(hd_obj) <- "Spatial.008um"
hd_obj <- RunBanksy(
  hd_obj,
  lambda = 0.2,
  verbose = FALSE,
  assay = "Spatial.008um",
  slot = "data",
  features = "variable",
  k_geom = 50
)
DefaultAssay(hd_obj) <- "BANKSY"
hd_obj <- RunPCA(hd_obj, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(hd_obj), npcs = 30)
hd_obj <- FindNeighbors(hd_obj, reduction = "pca.banksy", dims = 1:30)
hd_obj <- FindClusters(hd_obj, cluster.name = "banksy_cluster", resolution = 0.5)

saveRDS(hd_obj, "results/visiumhd_sketch_banksy.rds")
