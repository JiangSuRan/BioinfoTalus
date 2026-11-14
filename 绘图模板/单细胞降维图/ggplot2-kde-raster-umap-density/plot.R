## 暗背景等细胞下采样 UMAP 二维核密度分布图
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)

seurat_obj <- readRDS("Data/seurat_umap.rds")

coord <- Embeddings(seurat_obj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
meta <- cbind(seurat_obj@meta.data, coord)

# 跨条件按最小样本细胞数等量下采样，消除总细胞量差异导致的视觉密度偏差
min_n <- min(table(meta$condition))
set.seed(123)
meta_balanced <- meta %>%
  group_by(condition) %>%
  slice_sample(n = min_n) %>%
  ungroup()

p <- ggplot(meta_balanced, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  geom_point(color = "white", size = 0.05) +
  scale_fill_viridis(option = "magma") +
  facet_wrap(~ condition) +
  theme_grey(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = "black", color = NA),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey20", color = "grey10"),
    strip.text = element_text(size = 10, color = "white"),
    plot.background = element_rect(fill = "black", color = "black"),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

ggsave(filename = "results/umap_kde_density_dark.pdf", plot = p, width = 12, height = 5, units = "in")
