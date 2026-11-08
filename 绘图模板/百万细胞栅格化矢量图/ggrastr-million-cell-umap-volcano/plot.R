## ggrastr 百万细胞 UMAP 点层光栅化（矢量坐标轴 + 位图点层）
library(Seurat)
library(ggplot2)
library(ggrastr)
library(dplyr)

seu <- readRDS("Data/large_seurat.rds")

# 高密度点层 raster 化（raster.dpi 控制位图清晰度），坐标轴/文字保持矢量
p <- DimPlot(seu, reduction = "umap", group.by = "joint.state", raster = FALSE) +
  ggrastr::geom_point_rast(size = 0.05, raster.dpi = 1000) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 10)
  )

ggsave(filename = "results/ggrastr_umap_rasterized.pdf", plot = p, width = 8, height = 6, units = "in")
