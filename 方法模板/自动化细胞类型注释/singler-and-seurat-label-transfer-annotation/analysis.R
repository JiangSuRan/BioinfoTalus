## SingleR 参考集注释（含群内多数投票）与 Seurat 锚点标签映射（Label Transfer）
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(dplyr)

query_obj <- readRDS("Data/query_seurat.rds")
ref_obj <- readRDS("Data/reference_seurat.rds")
ref_se <- readRDS("Data/reference_singleR_se.rds")

# 方法 1：基于单细胞表达谱与参考集 Spearman 相关逐细胞推断并按 seurat_clusters 汇总主导标签，避免群均化丢失异质性
test_mat <- LayerData(query_obj, assay = "RNA", layer = "data")
singler_res <- SingleR(
  test = test_mat,
  ref = ref_se,
  labels = ref_se$label.main
)
query_obj$SingleR_cell_label <- singler_res$labels

cluster_vote <- data.frame(
  cluster = as.character(query_obj$seurat_clusters),
  label = as.character(singler_res$labels),
  stringsAsFactors = FALSE
) %>%
  count(cluster, label) %>%
  group_by(cluster) %>%
  slice_max(n, n = 1L, with_ties = FALSE) %>%
  ungroup()

vote_map <- setNames(cluster_vote$label, cluster_vote$cluster)
query_obj$SingleR_cluster_label <- unname(vote_map[as.character(query_obj$seurat_clusters)])

# 方法 2：基于 PCA 投影与共同锚点将自定义单细胞参考集标签与置信度逐细胞映射至查询集
anchors <- FindTransferAnchors(
  reference = ref_obj,
  query = query_obj,
  dims = 1:30,
  reference.reduction = "pca"
)
pred <- TransferData(anchorset = anchors, refdata = ref_obj$celltype, dims = 1:30)
query_obj <- AddMetaData(query_obj, metadata = pred)

saveRDS(query_obj, "results/annotated_query_seurat.rds")
write.csv(as.data.frame(singler_res), "results/singler_cell_scores.csv", row.names = TRUE)
write.csv(cluster_vote, "results/singler_cluster_votes.csv", row.names = FALSE)
