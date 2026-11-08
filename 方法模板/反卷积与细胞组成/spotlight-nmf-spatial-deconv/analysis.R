## SPOTlight seeded NMF 空间转录组反卷积与空间饼图权重导出
library(Seurat)
library(SPOTlight)

se_sc <- readRDS("Data/sc_reference_seurat.rds")
se_sp <- readRDS("Data/visium_seurat.rds")

Idents(se_sc) <- "subclass"
cluster_markers_all <- FindAllMarkers(se_sc, only.pos = TRUE)

# seeded NMF：单细胞参考的聚类标记作初始矩阵，反卷积得到 spot × 细胞类型权重
spotlight_decon <- spotlight_deconvolution(
  se_sc = se_sc,
  counts_spatial = se_sp@assays$Spatial@counts,
  clust_vr = "subclass",
  cluster_markers = cluster_markers_all,
  cl_n = 100,
  hvg = 3000,
  ntop = NULL,
  transf = "uv",
  method = "nsNMF",
  min_cont = 0.05
)

decon_mtrx <- spotlight_decon[[2]]
se_sp@meta.data <- cbind(se_sp@meta.data, decon_mtrx)

# 占比最高的类型作为 spot 的主导标签并写回
decon_mtrx <- as.data.frame(decon_mtrx)
decon_mtrx$spotlight.cluster <- colnames(decon_mtrx)[
  apply(decon_mtrx[, -ncol(decon_mtrx)], 1, which.max)
]
se_sp$spotlight.cluster <- factor(decon_mtrx$spotlight.cluster,
                                  levels = setdiff(colnames(decon_mtrx), "res_ss"))

write.csv(se_sp@meta.data, "results/spotlight_deconvolution_weights.csv")
saveRDS(se_sp, "results/spotlight_seurat.rds")

pdf("results/spotlight_scatterpie.pdf", width = 12, height = 10)
SPOTlight::spatial_scatterpie(
  se_obj = se_sp,
  cell_types_all = setdiff(colnames(decon_mtrx), "res_ss"),
  pie_scale = 0.4
) 
dev.off()
