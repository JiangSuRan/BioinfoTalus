## 多供体/样本分面 UMAP 亚群分布图（分面标注细胞总数）
library(dplyr)
library(ggplot2)

cell_df <- readRDS("Data/donor_umap_metadata.rds")

cell_df <- cell_df %>%
  add_count(donor, name = "n_cells") %>%
  mutate(donor_label = paste0(donor, "\n[", n_cells, "]"))

p <- ggplot(cell_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.8) +
  facet_wrap(~ donor_label, ncol = 4) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(filename = "results/donor_faceted_count_umap.pdf", plot = p, width = 10, height = 6, units = "in")
