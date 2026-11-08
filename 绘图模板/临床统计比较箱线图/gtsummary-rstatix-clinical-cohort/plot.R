## rstatix + ggpubr 多亚组精确检验标注的组间比较图
library(ggpubr)
library(rstatix)
library(ggplot2)
library(dplyr)

prop_df <- read.csv("Data/diff_abundance_props.csv", stringsAsFactors = FALSE)

# rstatix 精确组间 Wilcoxon + 位置计算
stat_test <- prop_df %>%
  wilcox_test(prop ~ progression_group) %>%
  add_xy_position(x = "celltype")

p <- ggboxplot(
  prop_df,
  x = "celltype",
  y = "prop",
  color = "progression_group",
  palette = "jco",
  outlier.shape = NA
) +
  stat_pvalue_manual(stat_test, label = "{scales::pvalue(p)}", size = 4, hide.ns = TRUE) +
  labs(x = NULL, y = "Cell type proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "results/rstatix_grouped_boxplot.pdf", plot = p, width = 8, height = 5, units = "in")
