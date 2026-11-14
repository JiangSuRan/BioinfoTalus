## CITE-seq 双模态 RNA + ADT 加权最近邻（WNN）整合与多模态 UMAP 聚类
library(Seurat)
library(dplyr)

bm <- readRDS("Data/citeseq_seurat.rds")

DefaultAssay(bm) <- "RNA"
n_pca <- min(30, ncol(bm) - 1)
bm <- NormalizeData(bm) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = n_pca, verbose = FALSE)

DefaultAssay(bm) <- "ADT"
VariableFeatures(bm) <- rownames(bm[["ADT"]])
n_apca <- min(18, nrow(bm[["ADT"]]) - 1, ncol(bm) - 1)
# 表面蛋白跨单细胞做 CLR 中心对数比标准化以消除细胞表面积与本底吸附差异
bm <- NormalizeData(bm, normalization.method = "CLR", margin = 2) %>%
  ScaleData() %>%
  RunPCA(reduction.name = "apca", npcs = n_apca, verbose = FALSE)

bm <- FindMultiModalNeighbors(
  bm,
  reduction.list = list("pca", "apca"),
  dims.list = list(seq_len(n_pca), seq_len(n_apca)),
  k.nn = min(20, ncol(bm) - 1),
  knn.range = min(200, ncol(bm) - 1),
  modality.weight.name = c("RNA.weight", "ADT.weight")
)

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

saveRDS(bm, "results/citeseq_wnn_seurat.rds")
