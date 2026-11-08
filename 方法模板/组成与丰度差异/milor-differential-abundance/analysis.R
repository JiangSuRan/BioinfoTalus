## miloR 单细胞 KNN 邻域构建与协变量校正的差异丰度检验
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(dplyr)

s.sub <- readRDS("Data/monocyte_seurat.rds")

# miloR 以 SingleCellExperiment 为容器；REF.SPCA 为上游已存储的监督降维结果
pbmc_small_sce <- as.SingleCellExperiment(s.sub)
traj_milo <- Milo(pbmc_small_sce)

traj_milo <- buildGraph(traj_milo, k = 10, d = 30, reduced.dim = "REF.SPCA")
# refined = TRUE 借助采样重定位使邻域均匀覆盖表型空间
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d = 30, reduced_dims = "REF.SPCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample = "sample")
traj_milo <- calcNhoodDistance(traj_milo, d = 30, reduced.dim = "REF.SPCA")

# 设计表为样本级协变量去重提取，行名须与 countCells 的样本名一致
traj_design <- data.frame(colData(traj_milo))[, c("sample", "class", "age")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample

# 年龄作为协变量纳入设计矩阵，避免 class 差异与年龄混杂
da_results <- testNhoods(traj_milo, reduced.dim = "REF.SPCA", design = ~ age + class, design.df = traj_design)

traj_milo <- buildNhoodGraph(traj_milo)
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")

saveRDS(traj_milo, "results/milo_object.rds")
write.csv(da_results, "results/milo_da_results.csv", row.names = FALSE)
