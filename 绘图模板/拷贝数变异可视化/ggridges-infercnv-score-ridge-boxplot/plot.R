## InferCNV 恶性 CNV 负荷亚群排序箱线图与尾概率梯度山脊图
library(dplyr)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(viridis)
library(patchwork)

cnv_mat <- readRDS("Data/infercnv_obs_mat.rds") %>% as.matrix()
cell_meta <- readRDS("Data/cell_meta.rds")

# 将各细胞观测值标准化并映射回 [-1, 1]，以偏移平方和量化全基因组拷贝数扰动总负荷
scaled_mat <- cnv_mat %>%
  t() %>%
  scale() %>%
  rescale(to = c(-1, 1)) %>%
  t()
cnv_scores <- colSums(scaled_mat * scaled_mat)

shared_cells <- intersect(names(cnv_scores), rownames(cell_meta))
plot_df <- data.frame(
  cnv_score = as.numeric(cnv_scores[shared_cells]),
  celltype = factor(cell_meta[shared_cells, "celltype"])
)

p_box <- ggplot(plot_df, aes(x = reorder(celltype, cnv_score, FUN = median), y = cnv_score, color = celltype)) +
  geom_boxplot(outlier.size = 0.7) +
  stat_compare_means() +
  labs(x = NULL, y = "CNV Scores", title = "CNV Burden by Cell Type") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11, angle = 45, vjust = 1, hjust = 1)
  )

p_ridge <- ggplot(plot_df, aes(x = cnv_score, y = reorder(celltype, cnv_score, FUN = median), fill = 0.5 - abs(0.5 - after_stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  labs(x = "CNV Score", y = NULL, title = "CNV Score ECDF Tail Density") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"))

p_all <- p_box / p_ridge
ggsave("results/infercnv_score_box_ridge.pdf", plot = p_all, width = 8, height = 8.5, units = "in")
