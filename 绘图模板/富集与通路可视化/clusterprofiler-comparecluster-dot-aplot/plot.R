## 多亚群 compareCluster 富集网络与顶部细胞数对齐气泡图
library(DOSE)
library(enrichplot)
library(ggplot2)
library(aplot)
library(dplyr)

kk <- readRDS("Data/compare_cluster_go.rds")
counts_df <- readRDS("Data/cluster_counts.rds")

kk_sim <- pairwise_termsim(kk)
p_emap <- emapplot(kk_sim, showCategory = 10)
ggsave("results/comparecluster_emapplot.pdf", plot = p_emap, width = 9, height = 8, units = "in")

p_dot <- dotplot(kk, label_format = 60) +
  scale_color_gradientn(colours = c("#E54924", "#498EA4")) +
  labs(x = "Cell Clusters", title = NULL) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9)
  )
# 剥离 dotplot 默认附加的 "(n=...)" 换行字样，使横轴分类与顶部细胞数柱图直接对齐
p_dot$data$Cluster <- sub("\n.*", "", as.character(p_dot$data$Cluster)) %>%
  factor(levels = as.character(counts_df$cluster))

p_top <- ggplot(counts_df, aes(x = factor(cluster, levels = cluster), y = n_cells, fill = factor(cluster, levels = cluster))) +
  geom_col(show.legend = FALSE) +
  labs(x = NULL, y = "Cell count") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 9)
  )

p_combo <- insert_top(p_dot, p_top, height = 0.22)
ggsave("results/comparecluster_dot_topcount.pdf", plot = p_combo, width = 8.5, height = 9.5, units = "in")
