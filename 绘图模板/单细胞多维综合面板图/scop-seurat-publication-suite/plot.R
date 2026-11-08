## scop 单细胞多维综合面板（CellDimPlot + GroupHeatmap 栅格化输出）
library(scop)
library(thisplot)
library(dplyr)

pancreas_sub <- readRDS("Data/scop_input.rds")
pancreas_sub <- standard_scop(pancreas_sub)

# 降维注释面板：分群着色 UMAP（带标签）
p1 <- CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

# 标志基因分组热图：按 CellType/SubCellType 双层聚合的均值表达
ht1 <- GroupHeatmap(
  pancreas_sub,
  features = c("Sox9", "Anxa2", "Bicc1", "Neurog3", "Hes6", "Fev", "Neurod1", "Rbp4", "Pyy", "Ins1", "Gcg", "Sst"),
  group.by = c("CellType", "SubCellType")
)

ggsave(filename = "results/scop_celldimplot.pdf", plot = p1, width = 6, height = 5, units = "in")

thisplot::panel_fix(
  ht1$plot,
  height = 4,
  width = 6,
  raster = TRUE,
  dpi = 300,
  filename = "results/scop_groupheatmap.pdf"
)
