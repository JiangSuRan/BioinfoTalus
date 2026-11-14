## 样本组排序免疫浸润相对比例堆叠柱状图与分组箱线比较联图
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(patchwork)

immune_mat <- readRDS("Data/ssgsea_immune_scores.rds")
sample_meta <- readRDS("Data/sample_group.rds")

# Min-Max 标准化后按样本归一化为相对构成，使堆叠柱反映免疫微环境结构比重
norm_mat <- apply(immune_mat, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-12)) %>%
  t()
rel_mat <- sweep(norm_mat, 2, colSums(norm_mat), FUN = "/")

sample_order_df <- sample_meta %>%
  filter(Sample %in% colnames(rel_mat)) %>%
  mutate(Group = factor(Group, levels = c("Control", "Treat"))) %>%
  arrange(Group, Sample)

stack_df <- rel_mat[, sample_order_df$Sample, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(ImmuneCell = rownames(.)) %>%
  melt(id.vars = "ImmuneCell", variable.name = "Sample", value.name = "Proportion") %>%
  left_join(sample_order_df, by = "Sample") %>%
  mutate(Sample = factor(Sample, levels = sample_order_df$Sample))

n_cells <- length(unique(stack_df$ImmuneCell))
cell_pal <- colorRampPalette(brewer.pal(12, "Paired"))(n_cells)

p_stack <- ggplot(stack_df, aes(x = Sample, y = Proportion, fill = ImmuneCell)) +
  geom_col(width = 0.88) +
  scale_fill_manual(values = cell_pal) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Relative Immune Proportion", fill = "Immune Cell") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 5, 0, 5)
  )

# 底部添加样本分组色带，清晰界定 Control 与 Treat 分组边界
p_bar <- ggplot(sample_order_df, aes(x = factor(Sample, levels = Sample), y = 1, fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Treat" = "#E64B35")) +
  labs(x = "Samples (ordered by Group)", y = "", fill = "Group") +
  theme_void(base_size = 10) +
  theme(
    axis.title.x = element_text(margin = margin(t = 4)),
    legend.position = "bottom",
    plot.margin = margin(0, 5, 5, 5)
  )

raw_df <- immune_mat[, sample_order_df$Sample, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(ImmuneCell = rownames(.)) %>%
  melt(id.vars = "ImmuneCell", variable.name = "Sample", value.name = "Score") %>%
  left_join(sample_order_df, by = "Sample")

p_box <- ggplot(raw_df, aes(x = ImmuneCell, y = Score, fill = Group)) +
  geom_boxplot(width = 0.65, outlier.size = 0.8, alpha = 0.85, position = position_dodge(0.75)) +
  stat_compare_means(aes(group = Group), method = "wilcox.test", label = "p.signif", hide.ns = FALSE) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Treat" = "#E64B35")) +
  labs(x = "", y = "ssGSEA Score", fill = "Group") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, colour = "black"),
    legend.position = "top"
  )

p_combined <- p_stack / p_bar / p_box +
  plot_layout(heights = c(10, 0.9, 9))

ggsave(filename = "results/immune_stacked_bar_group_boxplot.pdf", plot = p_combined, width = 11, height = 8.5, units = "in")
