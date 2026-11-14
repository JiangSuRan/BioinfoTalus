## 单细胞亚群组织偏好 Ro/e 观测期望比热图
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

meta <- readRDS("Data/cell_tissue_metadata.rds")

obs <- table(as.character(meta$cell_type), as.character(meta$tissue))
exp_mat <- outer(rowSums(obs), colSums(obs)) / sum(obs)
roe <- unclass(obs) / exp_mat

col_fun <- colorRamp2(c(0, 1, 1.5, 2), c("#fffde7", "#ffe0b2", "#ff9800", "#e65100"))

ht <- Heatmap(
  roe,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D",
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", roe[i, j]), x, y, gp = gpar(fontsize = 6))
  },
  width = ncol(roe) * unit(0.45, "inch"),
  height = nrow(roe) * unit(0.22, "inch"),
  name = "Ro/e"
)

pdf("results/roe_tissue_preference_heatmap.pdf", width = 4.5, height = 5.5)
draw(ht)
dev.off()
