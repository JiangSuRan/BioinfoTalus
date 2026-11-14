## SpaGene 空间可变基因检测、局部空间模式识别与配体-受体共定位推断
library(Seurat)
library(SpaGene)

st_obj <- readRDS("Data/visium_spatial.rds")
expr_mat <- LayerData(st_obj, assay = "Spatial", layer = "counts")
coords <- GetTissueCoordinates(st_obj)[, c("x", "y")]

# 基于最近邻空间图上的 EM 算法识别局部聚集高表达的空间模式基因（SVG）及多亚模式聚类
sg <- SpaGene(expr_mat, coords)
patterns <- FindPattern(sg, nPattern = 8)

# 结合配体-受体数据库筛选同一空间微环境邻域内显著共定位表达的细胞通讯分子对
data(LRpair_human)
lr_res <- SpaGene_LR(expr_mat, coords, LRpair = LRpair_human, knn = 8)

saveRDS(sg, "results/spagene_svg_results.rds")
saveRDS(patterns, "results/spagene_patterns.rds")
write.csv(lr_res, "results/spagene_colocalized_lr.csv", row.names = FALSE)
