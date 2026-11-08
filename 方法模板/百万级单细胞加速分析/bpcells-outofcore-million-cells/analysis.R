## BPCells 磁盘外存矩阵读取、Seurat v5 大图谱对象构建与归一化
library(BPCells)
library(Seurat)
library(dplyr)

# 打开 BPCells 目录格式矩阵（bit-packed 压缩，不进内存）
input_dir <- "Data/bpcells_hao2021"
count_matrix <- open_matrix_dir(dir = input_dir)

# Seurat v5 on-disk assay：counts 保留为 BPCells 引用，NormalizeData 支持流式
options(Seurat.object.assay.version = "v5")
seurat_object <- CreateSeuratObject(counts = count_matrix)

cell_metadata <- read.csv("Data/cell_metadata.csv")
rownames(cell_metadata) <- cell_metadata$CellID
seurat_object@meta.data <- cbind(seurat_object@meta.data, cell_metadata[colnames(seurat_object), ])

seurat_object <- NormalizeData(seurat_object)

saveRDS(seurat_object, "results/bpcells_seurat_v5.rds")

# 目录格式写出（供下游继续磁盘外存处理）
write_matrix_dir(
  mat = GetAssayData(seurat_object, assay = "RNA", layer = "counts"),
  dir = "results/bpcells_counts_out",
  overwrite = TRUE
)
