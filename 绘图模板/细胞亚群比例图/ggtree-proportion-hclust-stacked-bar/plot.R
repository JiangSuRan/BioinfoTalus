## 亚群构成层次聚类树（ggtree）与样本量标注堆叠比例对齐图
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtree)
library(patchwork)

meta <- readRDS("Data/cohort_subset_metadata.rds")

prop_df <- meta %>%
  group_by(cohort) %>%
  mutate(cohort_label = paste0(cohort, " (N = ", length(unique(sample_id)), ")")) %>%
  ungroup() %>%
  count(cohort_label, cell_type) %>%
  group_by(cohort_label) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

wide_mat <- prop_df %>%
  select(cohort_label, cell_type, pct) %>%
  pivot_wider(names_from = cell_type, values_from = pct, values_fill = 0) %>%
  as.data.frame()
rownames(wide_mat) <- wide_mat$cohort_label
wide_mat$cohort_label <- NULL

hc <- hclust(dist(wide_mat, method = "euclidean"), method = "average")
p_tree <- ggtree(hc) + geom_tiplab(size = 3) + hexpand(0.45)

# 直接从 ggtree 渲染坐标提取自底向上叶节点顺序，确保与 ladderize 后的聚类树严格逐行对齐
ordered_cohorts <- p_tree$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)
prop_df$cohort_label <- factor(prop_df$cohort_label, levels = ordered_cohorts)

p_bar <- ggplot(prop_df, aes(x = cohort_label, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  coord_flip() +
  labs(x = NULL, y = "Proportion", fill = "Subset") +
  theme_cowplot()

ggsave(filename = "results/hclust_subset_proportion_bar.pdf", plot = p_tree + p_bar + plot_layout(widths = c(1, 2)), width = 10, height = 5, units = "in")
