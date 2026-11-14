## SuperCell 图粗粒化 Metacell 构建、纯度评估与加权 PCA 转入 Seurat
library(SuperCell)
library(Seurat)
library(dplyr)
library(Matrix)
library(igraph)

sc_obj <- readRDS("Data/seurat_obj.rds")
counts_mat <- GetAssayData(sc_obj, assay = "RNA", layer = "counts")
cell_meta <- sc_obj$celltype

# gamma = 20 将样本规模压缩约 20 倍，大幅降低超大单细胞矩阵下游推断耗时
sc_meta <- SCimplify(counts_mat, k.knn = 5, gamma = 20, n.var.genes = 1000)
mc_ge <- supercell_GE(counts_mat, sc_meta$membership)
colnames(mc_ge) <- paste0("MC_", seq_len(ncol(mc_ge)))

sc_meta$cell_line <- supercell_assign(clusters = cell_meta, supercell_membership = sc_meta$membership, method = "jaccard")
sc_meta$purity <- supercell_purity(clusters = cell_meta, supercell_membership = sc_meta$membership, method = "entropy")

meta_df <- data.frame(
  size = sc_meta$supercell_size,
  cell_line = sc_meta$cell_line,
  purity = sc_meta$purity,
  row.names = colnames(mc_ge)
)

m_seu <- CreateSeuratObject(counts = mc_ge, meta.data = meta_df) %>%
  NormalizeData() %>%
  ScaleData(verbose = FALSE)

var_genes <- sort(sc_meta$genes.use)
VariableFeatures(m_seu) <- var_genes
n_comp <- min(30, ncol(mc_ge) - 1, length(var_genes) - 1)
m_seu <- RunPCA(m_seu, npcs = n_comp, verbose = FALSE)

# 必须按各 Metacell 实际覆盖的单细胞数目赋予 PCA 权重，而非等价权重
mc_pca <- supercell_prcomp(
  X = Matrix::t(mc_ge[var_genes, ]),
  genes.use = var_genes,
  fast.pca = TRUE,
  supercell_size = meta_df$size,
  k = ncol(Embeddings(m_seu, "pca"))
)
dimnames(mc_pca$x) <- dimnames(Embeddings(m_seu, "pca"))
m_seu[["pca"]]@cell.embeddings <- mc_pca$x
m_seu <- FindNeighbors(m_seu, dims = seq_len(n_comp), verbose = FALSE) %>%
  FindClusters(resolution = 0.8, verbose = FALSE) %>%
  RunUMAP(dims = seq_len(min(10, n_comp)), verbose = FALSE)

saveRDS(sc_meta, "results/supercell_obj.rds")
saveRDS(m_seu, "results/metacell_seurat.rds")
