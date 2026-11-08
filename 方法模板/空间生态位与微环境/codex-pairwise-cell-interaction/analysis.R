## CODEX 空间单细胞配对邻近互作富集统计
library(SpatialExperiment)
library(imcRtools)
library(data.table)

cur_features <- fread("Data/codex_features.csv", data.table = FALSE)
counts <- cur_features[, grepl("Cell: Mean", colnames(cur_features))]
meta <- cur_features[, c(1:5, 8:21, 762)]
coords <- cur_features[, c("Centroid X µm", "Centroid Y µm")]
colnames(coords) <- c("Pos_X", "Pos_Y")
cur_ch <- strsplit(colnames(counts), ":")
colnames(counts) <- sapply(cur_ch, function(x) x[1])

spe <- SpatialExperiment(
  assays = list(counts = t(counts)),
  colData = meta,
  sample_id = as.character(meta$Image),
  spatialCoords = as.matrix(coords)
)
colnames(spe) <- paste0(spe$sample_id, "_", 1:length(counts))

# asinh 变换将线性蛋白强度拉近似对称，缓解高表达通道主导距离与图结构
assay(spe, "exprs") <- asinh(counts(spe) / 1)
rowData(spe)$use_channel <- !grepl("DAPI", rownames(spe))

# 20 近邻图定义"邻近"关系，对每图像内细胞类型对做置换富集检验
spe <- buildSpatialGraph(spe, img_id = "Image", type = "knn", k = 20)
interaction_mat <- testInteractions(
  spe,
  group_by = "Image",
  label = "celltype",
  method = "classic",
  colPairName = "knn_interaction_graph",
  BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
)

write.csv(as.data.frame(interaction_mat), "results/codex_pairwise_interactions.csv", row.names = FALSE)
