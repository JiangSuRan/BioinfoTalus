## AUCell 通路活性 UMAP 质心标注与分组分面梯度图
library(Seurat)
library(AUCell)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(patchwork)

sc_obj <- readRDS("Data/seurat_obj.rds")
cells_auc <- readRDS("Data/cells_auc.rds")
target_set <- readRDS("Data/target_geneset.rds")[[1]]

umap_coords <- Embeddings(sc_obj, reduction = "umap") %>% as.data.frame()
colnames(umap_coords)[1:2] <- c("umap_1", "umap_2")

auc_vec <- getAUC(cells_auc)[target_set, colnames(sc_obj)]
plot_df <- sc_obj@meta.data %>%
  bind_cols(umap_coords) %>%
  mutate(AUC = as.numeric(auc_vec))

# 以各聚类坐标中位数而非均值计算标签质心，降低离散离群点对标注位置的拉扯
centroids <- plot_df %>%
  group_by(seurat_clusters) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p_global <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = AUC), size = 0.9) +
  scale_color_viridis(option = "E") +
  geom_label_repel(
    data = centroids,
    aes(label = seurat_clusters),
    size = 3.5,
    label.size = 0.4,
    segment.color = NA
  ) +
  labs(title = target_set, x = "UMAP 1", y = "UMAP 2") +
  theme_bw()

p_facet <- ggplot(plot_df, aes(x = umap_1, y = umap_2, color = AUC)) +
  geom_point(size = 0.8) +
  scale_color_viridis(option = "A") +
  facet_grid(. ~ group) +
  labs(title = paste(target_set, "by Group"), x = "UMAP 1", y = "UMAP 2") +
  theme_light(base_size = 11)

p_all <- p_global / p_facet + plot_layout(heights = c(1, 0.9))
ggsave("results/aucell_umap_centroid_facet.pdf", plot = p_all, width = 9, height = 8, units = "in")
