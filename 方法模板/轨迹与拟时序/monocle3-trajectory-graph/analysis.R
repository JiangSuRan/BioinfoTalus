## Monocle3 单细胞轨迹图学习（复用 Seurat UMAP 坐标）
library(Seurat)
library(monocle3)

seu_obj <- readRDS("Data/seurat_clustered.rds")

# monocle3 以原始计数构建 cell_data_set，基因元数据仅需基因名列
data <- GetAssayData(seu_obj, assay = "RNA", layer = "counts")
cell_metadata <- seu_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

# 用 Seurat 已有 UMAP 坐标替换 monocle3 自建嵌入，使轨迹图与既有聚类可视化保持一致
cds.embed <- SingleCellExperiment::reducedDim(cds, "UMAP")
int.embed <- Embeddings(seu_obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
SingleCellExperiment::reducedDim(cds, "UMAP") <- int.embed

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

saveRDS(cds, "results/monocle3_cds.rds")
