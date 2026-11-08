## Signac 单细胞染色质可及性 TF-IDF 归一化、SVD (LSI) 降维与 Harmony 整合聚类
library(Seurat)
library(Signac)
library(harmony)
library(dplyr)

counts_atac <- readRDS("Data/atac_peak_counts.rds")
cell_metadata <- read.csv("Data/cell_metadata.csv", row.names = 1, stringsAsFactors = FALSE)

chrom <- CreateChromatinAssay(counts = counts_atac)
seu_atac <- CreateSeuratObject(counts = chrom, assay = "ATAC", meta.data = cell_metadata)

seu_atac <- RunTFIDF(seu_atac)
seu_atac <- FindTopFeatures(seu_atac, min.cutoff = "q0")
seu_atac <- RunSVD(seu_atac)

# LSI 嵌入上做 Harmony 批次校正，无需重投影
lsi_coords <- Embeddings(seu_atac, reduction = "lsi")
lsi_harmony <- harmony::HarmonyMatrix(
  data_mat = lsi_coords,
  meta_data = as.data.frame(seu_atac[[]]),
  vars_use = "Batch",
  do_pca = FALSE
)
seu_atac[["harmony_lsi"]] <- CreateDimReducObject(
  embeddings = as.matrix(lsi_harmony),
  key = "harmonyLSI_",
  assay = "ATAC"
)

seu_atac <- FindNeighbors(seu_atac, reduction = "harmony_lsi", dims = 1:50)
seu_atac <- FindClusters(seu_atac, algorithm = 4, resolution = 0.5, random.seed = 1)

saveRDS(seu_atac, "results/signac_clustered.rds")
