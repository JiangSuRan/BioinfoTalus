## 空间转录组连续特征定制渐变色映射图（SpatialFeaturePlot + ggplot2）
library(Seurat)
library(patchwork)
library(purrr)

seu_merge <- readRDS("data/spatial_seurat.rds")

mg <- c("EGR1", "ENC1", "HOPX")

# scale.data 层配合 alpha 双值把低表达压为透明，突出连续表达梯度；坐标翻转与 y 反转适配大切片方向
tmp <- SpatialFeaturePlot(
  seu_merge,
  images = "slice1",
  stroke = 0,
  layer = "scale.data",
  pt.size.factor = 3,
  features = mg,
  alpha = c(0, 1),
  combine = FALSE
) %>%
  map(~. + coord_flip() + scale_y_reverse())

p <- wrap_plots(tmp, nrow = 1, byrow = TRUE)

ggsave(filename = "result/spatial_feature_panel.pdf", plot = p, width = 7, height = 6, units = "in")
