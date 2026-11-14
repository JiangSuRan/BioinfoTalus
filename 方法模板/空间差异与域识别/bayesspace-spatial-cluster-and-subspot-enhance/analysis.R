## BayesSpace Visium 空间域贝叶斯 MRF 聚类与亚 Spot 分辨率增强
library(SingleCellExperiment)
library(BayesSpace)

sce <- readVisium("Data/visium_spaceranger_out")
sce <- spatialPreprocess(sce, platform = "Visium", n.PCs = 15, n.HVGs = 2000, log.normalize = TRUE)
sce <- qTune(sce, qs = seq(2, 10), platform = "Visium", d = 15)

# 采用 t 分布误差马尔可夫随机场先验平滑邻域标签，抑制作图噪点并划定空间微结构域
sce <- spatialCluster(sce, q = 7, platform = "Visium", d = 15, init.method = "mclust", model = "t", nrep = 10000, burn.in = 100)
sce_enhanced <- spatialEnhance(sce, q = 7, platform = "Visium", d = 15, model = "t", nrep = 50000, burn.in = 100, save.chain = FALSE)

markers <- readLines("Data/spatial_marker_genes.txt")
markers <- intersect(markers, rownames(sce))
sce_enhanced <- enhanceFeatures(sce_enhanced, sce, feature_names = markers, nrounds = 0)

saveRDS(sce, "results/bayesspace_spot_clustered.rds")
saveRDS(sce_enhanced, "results/bayesspace_subspot_enhanced.rds")
