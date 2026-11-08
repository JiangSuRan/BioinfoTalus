## ggplot2 对角分割双三角矩阵点图（同时展示 eRegulon 转录因子表达量、靶基因集可及性/表达 AUC 与 RSS 特异性得分）
library(ggplot2)
library(ggnewscale)
library(dplyr)

rel_data <- read.csv("data/eregulon_activity.csv", stringsAsFactors = FALSE)

# 以格心 (x, y) 生成下/上三角两个三顶点多边形，NA 断开路径
split_tile <- function(x, y, fill1, fill2) {
  data.frame(
    x = c(x - 0.5, x + 0.5, x + 0.5, NA, x + 0.5, x - 0.5, x - 0.5, NA),
    y = c(y - 0.5, y - 0.5, y + 0.5, NA, y + 0.5, y + 0.5, y - 0.5, NA),
    fill = c(fill1, fill1, fill1, NA, fill2, fill2, fill2, NA),
    part = c(rep(1, 4), rep(2, 4))
  )
}

cell_levels <- levels(factor(rel_data$Cell_type))
regulon_levels <- levels(factor(rel_data$Regulon))
rel_data <- rel_data %>%
  mutate(
    Cell_num = as.numeric(factor(Cell_type, levels = cell_levels)),
    Regulon_num = as.numeric(factor(Regulon, levels = regulon_levels))
  )

# 下三角映射 AUC，上三角映射 TF 表达量（行内 z-score 后截断 ±2.5）
split_tiles <- do.call(rbind, lapply(seq_len(nrow(rel_data)), function(i) {
  split_tile(rel_data$Cell_num[i], rel_data$Regulon_num[i],
             rel_data$AUC_scale[i], rel_data$Expression_scale[i])
}))
split_tiles$Cell_type <- rel_data$Cell_type[rep(seq_len(nrow(rel_data)), each = 8)]
split_tiles$Regulon <- rel_data$Regulon[rep(seq_len(nrow(rel_data)), each = 8)]

p <- ggplot() +
  geom_polygon(
    data = subset(split_tiles, part == 1),
    aes(x = x, y = y, group = interaction(Cell_type, Regulon), fill = fill)
  ) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       limits = c(-2.5, 2.5), name = "AUC (z)") +
  ggnewscale::new_scale_fill() +
  geom_polygon(
    data = subset(split_tiles, part == 2),
    aes(x = x, y = y, group = interaction(Cell_type, Regulon), fill = fill)
  ) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "BrBG")),
                       limits = c(-2.5, 2.5), name = "TF expr (z)") +
  geom_point(data = rel_data, aes(x = Cell_num, y = Regulon_num, size = RSS),
             color = "black", shape = 23, fill = "white") +
  scale_radius(range = c(0.5, 3), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq_along(cell_levels), labels = cell_levels) +
  scale_y_continuous(breaks = seq_along(regulon_levels), labels = regulon_levels) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_equal()

ggsave(filename = "results/scenicplus_eregulon_split_tile.pdf", plot = p, width = 12, height = 8, units = "in")
