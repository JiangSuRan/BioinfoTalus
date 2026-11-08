## CHOIR 单细胞层级统计置换聚类与亚群识别
library(Seurat)
library(CHOIR)

seu <- readRDS("Data/seurat_input.rds")

alpha <- 0.05

# CHOIR 逐层做置换检验决定是否继续分裂亚群，Harmony 按 Batch 列校正批次后再统计检验
seu <- CHOIR(
  object = seu,
  alpha = alpha,
  batch_correction_method = "Harmony",
  batch_labels = "Batch",
  use_assay = "RNA",
  use_slot = "counts",
  n_cores = 4,
  random_seed = 1
)

# CHOIR 的 UMAP 需基于其层级树专用的距离度量计算，不能用通用降维替代
seu <- runCHOIRumap(seu)

clusters <- data.frame(
  CellID = rownames(seu@meta.data),
  Clusters = seu@meta.data[, paste0("CHOIR_clusters_", alpha)]
)
write.csv(clusters, "results/choir_clusters.csv", row.names = FALSE)

saveRDS(seu, "results/seurat_choir.rds")
