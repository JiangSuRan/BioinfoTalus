## 样本×亚群 100% 堆叠比例柱状图（RColorBrewer 定性扩展配色）
library(dplyr)
library(forcats)
library(RColorBrewer)
library(ggplot2)

meta <- readRDS("Data/sample_celltype_metadata.rds")

prop_df <- meta %>%
  mutate(
    sample = as_factor(sample),
    cluster = as_factor(cluster)
  ) %>%
  count(cluster, sample, name = "n")

qual_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_colors <- unlist(mapply(brewer.pal, qual_pals$maxcolors, rownames(qual_pals)))
cluster_colors <- setNames(qual_colors[seq_len(nlevels(prop_df$cluster))], levels(prop_df$cluster))

p <- ggplot(prop_df, aes(x = sample, y = n, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "", y = "Ratio", fill = "Cluster") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1, vjust = 1)
  )

ggsave(filename = "results/sample_cluster_stacked_proportion.pdf", plot = p, width = 6.5, height = 5, units = "in")
