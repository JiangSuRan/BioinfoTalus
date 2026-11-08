## 空间转录组细胞亚群原位空间分布与高亮映射图（SpatialDimPlot）
library(Seurat)
library(ggplot2)

seu <- readRDS("data/spatial_seurat.rds")

# 亚群原位映射：非高亮细胞统一灰化，凸显目标亚群的空间结构域分布
p <- SpatialDimPlot(
  seu,
  group.by = "domain",
  cells.highlight = colnames(seu)[seu$domain == "L4a"],
  cols.highlight = c("#e31a1c", "grey85"),
  pt.size.factor = 1.6,
  stroke = 0
) +
  theme(legend.position = "right")

ggsave(filename = "result/spatial_dim_panel.pdf", plot = p, width = 8, height = 6, units = "in")
