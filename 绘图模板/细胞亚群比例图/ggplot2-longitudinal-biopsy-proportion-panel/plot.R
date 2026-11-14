## 单例纵向多活检亚群堆叠比例、绝对细胞数折线与归一化分簇柱图组合面板
library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(scales)

meta <- readRDS("Data/patient_longitudinal_metadata.rds")
meta$timepoint <- as_factor(meta$timepoint)

count_df <- meta %>%
  count(timepoint, cell_type, name = "n_cells") %>%
  group_by(cell_type) %>%
  mutate(norm_frac = n_cells / sum(n_cells)) %>%
  ungroup()

p_stack <- ggplot(meta, aes(x = timepoint, fill = cell_type)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion", fill = "Cell type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_line <- ggplot(count_df, aes(x = timepoint, y = n_cells, color = cell_type, group = cell_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = NULL, y = "Number of cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_dodge <- ggplot(count_df, aes(x = timepoint, y = norm_frac, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = NULL, y = "Fraction across biopsies") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_total <- meta %>%
  count(timepoint, name = "total_cells") %>%
  ggplot(aes(x = timepoint, y = total_cells)) +
  geom_bar(stat = "identity", fill = "#5B7C99") +
  geom_text(aes(label = total_cells), vjust = -0.3, size = 3.2) +
  labs(x = NULL, y = "Total cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_panel <- ggarrange(p_stack, p_line, p_dodge, p_total, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
ggsave(filename = "results/longitudinal_biopsy_proportion_panel.pdf", plot = p_panel, width = 12, height = 8, units = "in")
