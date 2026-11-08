## Seurat + Harmony 多数据集整合与查询集锚点映射（TransferData / ProjectUMAP）
library(Seurat)
library(harmony)

seu_ref <- readRDS("Data/seurat_ref.rds")
seu_qry <- readRDS("Data/seurat_query.rds")
gene_used <- readLines("Data/hvg_genes.txt")

# 参考集在共享特征基因上降维，Harmony 消除 dataset 间批次效应后训练可投影的 UMAP 模型
seu_ref <- RunPCA(seu_ref, features = gene_used, npcs = 15, verbose = FALSE)
seu_ref <- RunHarmony(seu_ref, group.by.vars = "dataset")
seu_ref <- RunUMAP(seu_ref, reduction = "harmony", reduction.name = "harmony.umap", dims = 1:15, return.model = TRUE)

# pcaproject 锚点策略在大参考集上比 CCA 更稳健，锚点同时用于标签转移与嵌入整合
anchors <- FindTransferAnchors(reference = seu_ref, query = seu_qry, reduction = "pcaproject", reference.reduction = "pca", dims = 1:15, features = gene_used)

seu_qry <- TransferData(anchorset = anchors, reference = seu_ref, query = seu_qry, refdata = "meta.cluster", dims = 1:15)
seu_qry <- IntegrateEmbeddings(anchorset = anchors, reference = seu_ref, query = seu_qry, new.reduction.name = "ref.pca", dims.to.integrate = 1:15)
seu_qry <- ProjectUMAP(query = seu_qry, query.reduction = "ref.pca", query.dims = 1:15, reference = seu_ref, reference.reduction = "pca", reference.dims = 1:15, reduction.model = "harmony.umap")

saveRDS(seu_qry, "results/seurat_query_mapped.rds")
