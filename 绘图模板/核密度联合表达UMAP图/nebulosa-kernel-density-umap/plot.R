## Nebulosa 双标志物联合核密度 UMAP 图（dropout 遮挡的可视化修正）
library(Seurat)
library(Nebulosa)

seu <- readRDS("Data/neutrophil_seurat.rds")
DimPlot(seu, reduction = "tsne", group.by = "celltype")

# 三算法中性粒细胞打分的联合核密度：核密度平滑解决稀疏共表达点被单点遮挡
p1 <- Nebulosa::plot_density(seu, "NEU_xcell", "NEU_mcp", "NEU_quantiseq",
                             size = 0.5, reduction = "tsne",
                             joint = TRUE) &
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p2 <- Nebulosa::plot_density(seu, "NEU_consensus", size = 0.5, reduction = "tsne") &
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank(), panel.grid.major = element_blank())

ggsave(filename = "results/nebulosa_joint_density.pdf", plot = p1, width = 10, height = 4, units = "in")
ggsave(filename = "results/nebulosa_single_density.pdf", plot = p2, width = 6, height = 5, units = "in")
