## imcRtools/cytomapper CODEX 空间成像邻接图构建、互作置换检验与细胞社区聚合
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
colnames(spe) <- paste0(spe$sample_id, "_", seq_len(length(counts)))
assay(spe, "exprs") <- asinh(counts(spe) / 1)
rowData(spe)$use_channel <- !grepl("DAPI", rownames(spe))

# KNN 邻接图 + 邻域互作置换检验
spe <- buildSpatialGraph(spe, img_id = "Image", type = "knn", k = 20)
interaction_out <- testInteractions(
  spe,
  group_by = "Image",
  label = "celltype",
  method = "classic",
  colPairName = "knn_interaction_graph",
  BPPARAM = BiocParallel::SerialParam(RNGseed = 123)
)

# 细胞社区：按邻域标签构成聚合检测共同出现的细胞类型组合
spe <- aggregateNeighbors(
  spe,
  colPairName = "knn_interaction_graph",
  aggregate_by = "metadata",
  count_by = "celltype",
  name = "aggregatedNeighbors"
)

community_out <- as.data.frame(colData(spe)[, grep("^aggregatedNeighbors", colnames(colData(spe)))])
community_out$Image <- spe$Image
community_out$celltype <- spe$celltype
write.csv(as.data.frame(interaction_out), "results/codex_interactions.csv", row.names = FALSE)
write.csv(community_out, "results/codex_neighbor_aggregation.csv", row.names = FALSE)
