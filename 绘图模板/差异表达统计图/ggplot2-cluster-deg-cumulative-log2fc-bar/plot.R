## 单细胞多亚群差异基因双向累计 log2FC 与上/下调数目统计条形图
library(dplyr)
library(forcats)
library(ggplot2)

deg_df <- readRDS("Data/cluster_markers.rds")

# 按亚群统计显著上下调基因的累计 log2FC 负荷，以净效应方向与强度排列亚群
deg_sig <- deg_df %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  mutate(Direction = if_else(avg_log2FC > 0, "Up", "Down"))

sum_df <- deg_sig %>%
  group_by(cluster, Direction) %>%
  summarise(
    cum_log2FC = sum(avg_log2FC),
    n_genes = n(),
    .groups = "drop"
  )

order_levels <- sum_df %>%
  group_by(cluster) %>%
  summarise(net_fc = sum(cum_log2FC), .groups = "drop") %>%
  arrange(net_fc) %>%
  pull(cluster)

sum_df <- sum_df %>%
  mutate(
    cluster = factor(cluster, levels = order_levels),
    label_y = if_else(cum_log2FC >= 0, cum_log2FC + max(abs(sum_df$cum_log2FC)) * 0.04, cum_log2FC - max(abs(sum_df$cum_log2FC)) * 0.04)
  )

p <- ggplot(sum_df, aes(x = cluster, y = cum_log2FC, fill = Direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  geom_text(aes(y = label_y, label = n_genes), size = 3.1) +
  scale_fill_manual(values = c("Up" = "#D6604D", "Down" = "#4393C3")) +
  coord_flip() +
  labs(x = "Cell Cluster", y = "Cumulative avg_log2FC (n DEGs labeled)", fill = "Regulation") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"))

ggsave(filename = "results/cluster_deg_cumulative_log2fc.pdf", plot = p, width = 7.5, height = 5.5, units = "in")
