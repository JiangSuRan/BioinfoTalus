## 生物学重复样本亚群占比中位数+IQR误差棒及散点显著性分面柱状图
library(dplyr)
library(ggplot2)

meta <- readRDS("Data/sample_treatment_metadata.rds")
stats_df <- readRDS("Data/celltype_treatment_fdr.rds")

prop_df <- meta %>%
  count(sample_id, treatment, cell_type) %>%
  group_by(sample_id) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

summary_df <- prop_df %>%
  group_by(treatment, cell_type) %>%
  summarise(
    median = median(pct),
    q25 = quantile(pct, 0.25),
    q75 = quantile(pct, 0.75),
    .groups = "drop"
  ) %>%
  left_join(stats_df, by = c("treatment", "cell_type")) %>%
  mutate(
    signif_label = case_when(
      fdr < 0.0001 ~ "****",
      fdr < 0.001 ~ "***",
      fdr < 0.01 ~ "**",
      fdr < 0.1 ~ "*",
      TRUE ~ ""
    )
  )

p <- ggplot(summary_df, aes(x = treatment, y = median, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2, position = position_dodge(0.9)) +
  geom_jitter(data = prop_df, aes(x = treatment, y = pct), pch = 21, alpha = 0.6, width = 0.2, size = 1.5) +
  geom_text(aes(x = treatment, y = q75 + max(prop_df$pct) * 0.05, label = signif_label), size = 5, fontface = "bold") +
  facet_grid(~ cell_type) +
  labs(x = NULL, y = "% of Sequenced Cells") +
  theme(
    axis.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = "grey94"),
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave(filename = "results/replicate_proportion_median_iqr.pdf", plot = p, width = 12, height = 4.5, units = "in")
