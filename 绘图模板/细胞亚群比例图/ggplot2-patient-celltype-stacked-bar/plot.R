## 按临床分组排序的患者级亚群堆叠比例柱状图（患者标签按分组着色）
library(dplyr)
library(ggplot2)
library(scales)
library(ggtext)

meta <- readRDS("Data/patient_celltype_metadata.rds")
condition_colors <- readRDS("Data/condition_palette.rds")

prop_df <- meta %>%
  mutate(
    condition = factor(condition, levels = names(condition_colors)),
    label_md = paste0("<span style='color:", condition_colors[as.character(condition)], "'>", patient_id, "</span>")
  ) %>%
  arrange(condition, patient_id) %>%
  mutate(label_md = factor(label_md, levels = unique(label_md))) %>%
  count(label_md, cell_type)

p <- ggplot(prop_df, aes(x = label_md, y = n, fill = cell_type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = "Patient", y = "Cell proportion", fill = "Cell type") +
  theme_classic() +
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1, size = 7)
  )

ggsave(filename = "results/patient_celltype_stacked_bar.pdf", plot = p, width = 12, height = 5, units = "in")
