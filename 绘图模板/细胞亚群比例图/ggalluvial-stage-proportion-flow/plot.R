## 发育阶段/时间进程亚群比例冲击流动堆积图（ggalluvial）
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggalluvial)
library(scales)

meta <- readRDS("Data/cell_stage_metadata.rds")
meta$stage <- as_factor(meta$stage)

prop_df <- meta %>%
  count(stage, cell_type) %>%
  complete(stage, cell_type, fill = list(n = 0)) %>%
  group_by(stage) %>%
  mutate(rate = n / sum(n)) %>%
  ungroup()

p <- ggplot(
  prop_df,
  aes(x = stage, stratum = cell_type, alluvium = cell_type, y = rate, fill = cell_type)
) +
  geom_flow() +
  geom_stratum(width = 0.6, colour = "grey30") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell proportion", fill = "Subtype") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.text = element_text(size = 11)
  )

ggsave(filename = "results/stage_proportion_alluvial.pdf", plot = p, width = 8, height = 6, units = "in")
