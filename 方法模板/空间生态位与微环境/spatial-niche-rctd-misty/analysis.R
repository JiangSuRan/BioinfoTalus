## 空间转录组细胞丰度矩阵构建与空间生态位（Niche）聚类表征
library(Seurat)
library(dplyr)

st_obj <- readRDS("Data/visium_seurat.rds")
rctd <- readRDS("Data/rctd_visium.rds")

# RCTD 权重（spot × 细胞类型）即每个 spot 的微环境组成，先与空间对象对齐并写回 metadata
weights.df <- as.data.frame(as.matrix(rctd@results$weights))
commonspots <- intersect(colnames(st_obj), rownames(weights.df))
st_obj <- subset(st_obj, cells = commonspots)
st_obj <- AddMetaData(st_obj, metadata = weights.df)

# 以细胞组成权重为特征的 Niche assay：归一化缩放后聚类，簇即为微环境生态位
st_obj[["niche"]] <- CreateAssayObject(counts = t(as.matrix(weights.df[, colnames(st_obj), drop = FALSE])))
DefaultAssay(st_obj) <- "niche"
st_obj <- NormalizeData(st_obj)
st_obj <- ScaleData(st_obj, features = rownames(st_obj))
st_obj <- RunPCA(st_obj, reduction.name = "niche.pca", verbose = FALSE)
st_obj <- FindNeighbors(st_obj, reduction = "niche.pca", dims = 1:10)
st_obj <- FindClusters(st_obj, resolution = 0.4)
st_obj$niche <- Idents(st_obj)

saveRDS(st_obj, "results/spatial_niche_seurat.rds")
