## 图谱队列细胞数平方根柱图、组织来源堆叠柱图与环状图组合面板
library(dplyr)
library(ggplot2)
library(cowplot)
library(scales)
library(patchwork)

meta <- readRDS("Data/cohort_cell_metadata.rds")

count_df <- meta %>%
  count(cohort, name = "cell_count")

p_count <- ggplot(count_df, aes(x = cohort, y = cell_count)) +
  geom_bar(stat = "identity", fill = "#93B5C6") +
  scale_y_continuous(trans = "sqrt", labels = comma) +
  labs(x = NULL, y = "Number of cells") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

source_df <- meta %>%
  count(cohort, tissue)

p_source <- ggplot(source_df, aes(x = cohort, y = n, fill = tissue)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = NULL, y = "Proportion", fill = "Tissue") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

donut_df <- meta %>%
  count(cohort, name = "count") %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2,
    label = paste0(cohort, " (", count, ")")
  )

p_donut <- ggplot(donut_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = cohort)) +
  geom_rect(color = "white", linewidth = 0.2) +
  geom_text(x = 3.5, aes(y = label_pos, label = label), size = 2.1) +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

p_all <- p_count + p_source + p_donut + plot_layout(ncol = 3, widths = c(1.2, 1.2, 1.1))
ggsave(filename = "results/cohort_composition_overview_panel.pdf", plot = p_all, width = 12.5, height = 4, units = "in")
