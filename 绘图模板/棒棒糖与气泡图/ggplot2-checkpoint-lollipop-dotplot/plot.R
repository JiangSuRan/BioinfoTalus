## 免疫检查点配受体表达棒棒糖图（ggplot2 双向/分级 Lollipop）
library(ggpubr)
library(dplyr)
library(insight)

all_data <- read.csv("data/checkpoint_dge.csv", stringsAsFactors = FALSE)

# 显著性星号标记，校正后 p>0.05 记 ns 并在图中剔除
all_data$sign <- insight::format_p(all_data$p_val_adj, stars_only = TRUE)
all_data <- all_data %>% mutate(sign = ifelse(p_val_adj > 0.05, "ns", sign))

plot_data <- all_data %>%
  mutate(gene = as.factor(gene), comparison = as.factor(comparison)) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(sign != "ns")

map_colours <- c("#b2182b", "#1F618D", "#f4a582", "#9B59B6", "#85929E", "#1c9099", "#74add1", "#053061",
                 "#1b7837", "#b8e186", "#bebada", "#fed976", "#e7298a", "#47E3FF", "#FFBF47",
                 "#A93226", "#270e26", "orange", "#b8bc53", "#5628ce", "#fa909c",
                 "#8ff331", "#FF6347", "#6347FF", "#556270", "#4ECDC4", "#C7F464", "#FF6B6B",
                 "#C44D58", "#E3FF47", "#FF4787")

# ggdotchart 以 segments + 大圆点构成棒棒糖图层，圆点内嵌显著性星号
p <- ggdotchart(
  plot_data,
  x = "gene",
  y = "avg_log2FC",
  color = "celltype",
  palette = map_colours,
  sorting = "descending",
  add = "segments",
  add.params = list(color = "lightgray", size = 1),
  group = "comparison",
  dot.size = 6,
  label = "sign",
  font.label = list(color = "white", size = 7, vjust = 0.5),
  ggtheme = theme_pubr()
) +
  theme_cleveland() +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  coord_cartesian(ylim = c(-1.4, 1.4)) +
  facet_wrap(~comparison, nrow = 1, scales = "free_x") +
  font("xy.text", size = 9, color = "black", face = "bold.italic")

ggsave(filename = "result/checkpoint_lollipop.pdf", plot = p, width = 8, height = 6, units = "in")
