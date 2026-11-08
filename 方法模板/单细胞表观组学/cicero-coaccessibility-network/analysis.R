## Cicero 聚合单细胞 ATAC 顺式共可及性网络（LSI/UMAP + run_cicero）
# 输入为 monocle3 cell_data_set（或由 Seurat counts 经 new_cell_data_set 构建）
library(cicero)
library(monocle3)

set.seed(2017)

input_cds <- readRDS("Data/atac_cds.rds")

input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP", preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

# 以 UMAP 邻域聚合细胞形成 kNN “伪细胞”矩阵，扩大共可及性估计的功效
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# 染色体长度表（mm10；人源换 hg19/hg38 对应文件）
chromosome_length <- read.table("Data/mm10_chromosome_length.txt")

# 图 Lasso 距离惩罚估计顺式共可及性连接（启动子-增强子配对的基础）
conns <- run_cicero(cicero_cds, chromosome_length)

write.csv(conns, "results/cicero_connections.csv", row.names = FALSE)
saveRDS(conns, "results/cicero_conns.rds")
