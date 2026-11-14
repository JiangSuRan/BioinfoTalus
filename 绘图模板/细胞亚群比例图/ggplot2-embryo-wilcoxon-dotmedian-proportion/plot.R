## 样本级亚群比例散点、中位数横线与成对 Wilcoxon 检验显著性图
library(dplyr)
library(forcats)
library(ggplot2)
library(ggsignif)

meta <- readRDS("Data/sample_period_metadata.rds")
meta$period <- as_factor(meta$period)

ratio_df <- meta %>%
  count(sample_id, period, cell_type) %>%
  group_by(sample_id, period) %>%
  mutate(ratio = n / sum(n)) %>%
  ungroup()

median_df <- ratio_df %>%
  group_by(cell_type, period) %>%
  summarise(median_ratio = median(ratio), .groups = "drop") %>%
  mutate(
    x_start = as.numeric(period) - 0.3,
    x_end = as.numeric(period) + 0.3
  )

comparisons <- combn(levels(meta$period), 2, simplify = FALSE)

p <- ggplot(ratio_df, aes(x = period, y = ratio, color = period)) +
  geom_jitter(width = 0.15, size = 1.8) +
  geom_segment(
    data = median_df,
    aes(x = x_start, xend = x_end, y = median_ratio, yend = median_ratio),
    inherit.aes = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.1,
    color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(x = NULL, y = "Cell proportion per sample") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave(filename = "results/sample_proportion_wilcoxon_dotmedian.pdf", plot = p, width = 10, height = 6, units = "in")
