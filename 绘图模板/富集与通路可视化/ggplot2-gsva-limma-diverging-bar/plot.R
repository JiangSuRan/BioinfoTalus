## GSVA 差异通路 limma t 值双向发散条形图
library(dplyr)
library(stringr)
library(ggplot2)
library(ggprism)

diff_tbl <- readRDS("Data/gsva_limma_diff.rds")

dat_plot <- data.frame(id = rownames(diff_tbl), t = diff_tbl$t) %>%
  mutate(
    id = str_replace(id, "^HALLMARK_", ""),
    threshold = factor(
      ifelse(t >= 1, "Up", ifelse(t <= -1, "Down", "NoSignifi")),
      levels = c("Up", "Down", "NoSignifi")
    )
  ) %>%
  arrange(t) %>%
  mutate(id = factor(id, levels = id))

# 条件子集替代下标区间切片，避免某一阈区间样本数为 0 时触发逆序索引错误
neg_df <- dat_plot %>% filter(t < 0) %>% mutate(txt_col = ifelse(t < -1, "black", "grey45"))
pos_df <- dat_plot %>% filter(t >= 0) %>% mutate(txt_col = ifelse(t >= 1, "black", "grey45"))

p <- ggplot(dat_plot, aes(x = id, y = t, fill = threshold)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#36638a", "NoSignifi" = "#cccccc", "Down" = "#7bcd7b")) +
  geom_hline(yintercept = c(-2, 2), color = "white", linewidth = 0.5, linetype = "dashed") +
  geom_text(data = neg_df, aes(x = id, y = 0.1, label = id, color = txt_col), hjust = 0, size = 3, show.legend = FALSE) +
  geom_text(data = pos_df, aes(x = id, y = -0.1, label = id, color = txt_col), hjust = 1, size = 3, show.legend = FALSE) +
  scale_color_identity() +
  labs(x = NULL, y = "t value of GSVA score") +
  guides(fill = "none") +
  theme_prism(border = TRUE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave("results/gsva_diverging_bar.pdf", plot = p, width = 8.5, height = max(5, 0.25 * nrow(dat_plot) + 1.5), units = "in")
