## 生物学重复亚群占比均值±标准误误差棒与样本散点分面柱状图
library(dplyr)
library(forcats)
library(ggplot2)

meta <- readRDS("Data/replicate_celltype_metadata.rds") %>%
  mutate(
    sample = as_factor(sample),
    group = as_factor(group),
    cell_type = as_factor(cell_type)
  )

# 先用完整交叉列联表算样本内各亚群比例，保证某重复中缺失的亚群计为 0 而非被剔除
sample_meta <- meta %>%
  distinct(sample, group)

tab <- table(meta$sample, meta$cell_type)
prop_mat <- prop.table(tab, margin = 1)
prop_df <- as.data.frame(as.table(prop_mat), responseName = "prop") %>%
  rename(sample = Var1, cell_type = Var2) %>%
  left_join(sample_meta, by = "sample")

p <- ggplot(prop_df, aes(x = group, y = prop, fill = group)) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(0.9)) +
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,
    position = position_dodge(0.9),
    width = 0.2
  ) +
  geom_point(size = 2.5, pch = 19) +
  facet_grid(~ cell_type) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Group", y = "Proportion", fill = "Group") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black")
  )

ggsave(filename = "results/replicate_mean_se_proportion_bar.pdf", plot = p, width = 10, height = 4, units = "in")
