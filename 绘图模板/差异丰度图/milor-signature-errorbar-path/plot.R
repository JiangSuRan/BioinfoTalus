## 多条件 Milo 差异丰度特征均值±标准差折线误差棒图
library(dplyr)
library(ggplot2)

da_table <- readRDS("Data/milo_multigroup_da.rds")

da_summary <- da_table %>%
  filter(!is.na(cell_type)) %>%
  group_by(contrast, cell_type) %>%
  summarise(mean = mean(logFC), sd = sd(logFC), .groups = "drop") %>%
  mutate(cell_type = factor(cell_type, levels = rev(unique(da_table$cell_type))))

p <- ggplot(da_summary, aes(y = cell_type, x = mean, color = contrast, group = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), orientation = "y", linewidth = 0.4, width = 0.2, alpha = 0.75) +
  geom_path(aes(group = contrast), linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(limits = c(-5, 5)) +
  labs(x = "Mean logFC (Milo signature)", y = NULL, color = "Contrast") +
  theme_minimal()

ggsave(filename = "results/milo_signature_errorbar_path.pdf", plot = p, width = 8, height = 10, units = "in")
