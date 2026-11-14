## Giotto 空间邻近网络构建、BinSpect 空间可变基因检测与 HMRF 空间域分割
library(Giotto)

g_obj <- readRDS("Data/giotto_object.rds")

g_obj <- normalizeGiotto(gobject = g_obj, scalefactor = 6000, verbose = FALSE)
g_obj <- calculateHVG(gobject = g_obj, show_plot = FALSE)
gene_meta <- fDataDT(g_obj)
feat_genes <- gene_meta[hvg == "yes"]$feat_ID

g_obj <- runPCA(gobject = g_obj, genes_to_use = feat_genes, scale_unit = FALSE)
g_obj <- createNearestNetwork(gobject = g_obj, dimensions_to_use = 1:10, k = 15)
g_obj <- doLeidenCluster(gobject = g_obj, resolution = 0.4, n_iterations = 1000)

# 构建 Delaunay 空间物理网络并用 BinSpect-kmeans 检测具有显著空间邻域相关性基因
g_obj <- createSpatialNetwork(gobject = g_obj, method = "Delaunay", name = "Delaunay_network")
km_svg <- binSpect(g_obj, calc_hub = TRUE, hub_min_int = 5, spatial_network_name = "Delaunay_network")
ext_cor <- detectSpatialCorFeats(g_obj, method = "network", spatial_network_name = "Delaunay_network", subset_feats = km_svg$feats[1:300])

# HMRF 基于空间关联基因对与局部邻域先验联合划分连续组织微环境空间域
hmrf_res <- doHMRF(gobject = g_obj, expression_values = "scaled", spatial_network_name = "Delaunay_network", spatial_genes = km_svg$feats[1:100], k = 6, betas = c(0, 10, 3))
g_obj <- addHMRF(gobject = g_obj, HMRFoutput = hmrf_res, k = 6, betas_to_add = c(0, 10, 20))

saveRDS(g_obj, "results/giotto_hmrf_analyzed.rds")
write.csv(km_svg, "results/giotto_binspect_svgs.csv", row.names = FALSE)
