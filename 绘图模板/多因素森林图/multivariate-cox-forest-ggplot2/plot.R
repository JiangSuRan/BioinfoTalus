## 多因素森林图（NC 风格多变量 Cox 回归，ggplot2 像素坐标手排布局）
library(ggplot2)

estimates <- read.csv(
  "data/nc_forest_estimates.csv",
  stringsAsFactors = FALSE,
  colClasses = c(p_label = "character", stars = "character")
)
cell_types <- read.csv("data/nc_cell_types.csv", stringsAsFactors = FALSE)

estimates  <- estimates[order(estimates$row_order), ]
cell_types <- cell_types[order(cell_types$legend_order), ]

estimates$stars[is.na(estimates$stars)] <- ""

cluster_colors <- c(
  "1"  = "#388B4F", "2"  = "#7FB864", "3"  = "#F5CCCC",
  "4"  = "#D7589C", "5"  = "#A59ACB", "6"  = "#6F649B",
  "7"  = "#006DDC", "8"  = "#6EB6FF", "9"  = "#B6D9FA",
  "10" = "#E0A7CB", "11" = "#8C4601"
)

forest_x <- function(x) {
  218 + (log(x) - log(0.5)) / (log(2.5) - log(0.5)) * (397 - 218)
}

estimates$y <- ifelse(
  estimates$row_order <= 3,
  292 - (estimates$row_order - 1) * 19,
  216 - (estimates$row_order - 4) * 19
)

estimates$x    <- forest_x(estimates$hr)
estimates$xmin <- forest_x(estimates$ci_low)
estimates$xmax <- forest_x(estimates$ci_high)

clinical   <- estimates[is.na(estimates$cluster_id), ]
signatures <- estimates[!is.na(estimates$cluster_id), ]

shade <- ifelse(
  estimates$row_order <= 3,
  estimates$row_order %% 2 == 1,
  (estimates$row_order - 3) %% 2 == 1
)

group_labels <- data.frame(
  label = c("Mono", "TAM", "DC"),
  x = 14,
  y = c(206.5, 140, 64)
)
group_lines <- data.frame(
  ymin = c(188.5, 95.5, 35.5),
  ymax = c(225.5, 184.5, 92.5)
)

cell_types$x <- c(505, 505, 505, 642, 505, 642, 505, 505, 505, 654, 505)
cell_types$y <- c(297, 278, 236, 236, 217, 217, 198, 157, 116, 116,  97)

cell_types$plot_label <- c(
  "plain('Classical monocyte')",
  "plain('Nonclassical monocyte')",
  "plain('FOLR2')^'+'~plain('TAM')",
  "plain('IFN-TAM')",
  "plain('Inflam-TAM')",
  "plain('RTM-TAM')",
  "plain('TREM2')^'+'~plain('TAM')",
  "plain('Cycling cell')",
  "plain('cDC1 (CLEC9A')^'+'*plain(')')",
  "plain('cDC2 (CD1C')^'+'*plain(')')",
  "plain('cDC3 (LAMP3')^'+'*plain(')')"
)

legend_headers <- data.frame(
  label = c("Monocyte:", "TAM:", "Proliferating:", "DC:"),
  x = 498,
  y = c(318, 257, 178, 137)
)

axis_data <- data.frame(
  label = c("0.5", "1.0", "1.5", "2.0", "2.5"),
  x = forest_x(c(0.5, 1, 1.5, 2, 2.5))
)

p <- ggplot() +

  geom_rect(
    data = estimates[shade, ],
    aes(xmin = 29, xmax = 491, ymin = y - 9.5, ymax = y + 9.5),
    fill = "#EEF2F1",
    colour = NA
  ) +

  geom_segment(
    data = group_lines,
    aes(x = 28, xend = 28, y = ymin, yend = ymax),
    linewidth = 0.35,
    colour = "#4A4A4A"
  ) +
  geom_text(
    data = group_labels,
    aes(x = x, y = y, label = label),
    angle = 90,
    size = 2.85
  ) +

  annotate(
    "segment", x = forest_x(1), xend = forest_x(1), y = 35, yend = 309,
    linewidth = 0.34, linetype = "22", colour = "#4B4B4B"
  ) +
  geom_segment(
    data = estimates,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 0.38,
    colour = "#333333"
  ) +
  geom_point(
    data = estimates,
    aes(x = x, y = y),
    shape = 15,
    size = 2.45
  ) +

  annotate(
    "text", x = 36, y = 315, label = "Variable",
    hjust = 0, fontface = "bold", size = 3.05
  ) +
  annotate(
    "text", x = 111, y = 315, label = "HR (95% CI)",
    hjust = 0, fontface = "bold", size = 3.05
  ) +
  annotate(
    "text", x = 415, y = 315, label = "bolditalic(P)~bold(value)",
    parse = TRUE, hjust = 0, size = 3.05
  ) +

  geom_text(
    data = clinical,
    aes(x = 36, y = y, label = variable),
    hjust = 0,
    size = 2.85
  ) +
  annotate(
    "text", x = 36, y = 235, label = "Cell signature",
    hjust = 0, size = 2.85
  ) +

  geom_point(
    data = signatures,
    aes(x = 51, y = y, fill = factor(cluster_id)),
    shape = 21, size = 4.2, stroke = 0, colour = "white"
  ) +
  geom_text(
    data = signatures,
    aes(x = 51, y = y, label = cluster_id),
    colour = "white", fontface = "bold", size = 2.05
  ) +

  geom_text(
    data = estimates,
    aes(x = 111, y = y, label = hr_ci_label),
    hjust = 0, size = 2.75
  ) +
  geom_text(
    data = estimates,
    aes(x = 415, y = y, label = p_label),
    hjust = 0, fontface = "italic", size = 2.75
  ) +
  geom_text(
    data = estimates[nzchar(estimates$stars), ],
    aes(x = 462, y = y, label = stars),
    hjust = 0, fontface = "bold", size = 2.55
  ) +

  annotate(
    "segment", x = 218, xend = 397, y = 30, yend = 30,
    linewidth = 0.35, colour = "#4A4A4A"
  ) +
  geom_segment(
    data = axis_data,
    aes(x = x, xend = x, y = 30, yend = 26),
    linewidth = 0.35, colour = "#4A4A4A"
  ) +
  geom_text(
    data = axis_data,
    aes(x = x, y = 15, label = label),
    colour = "#333333", size = 2.65
  ) +

  geom_text(
    data = legend_headers,
    aes(x = x, y = y, label = label),
    hjust = 0, fontface = "bold", size = 3
  ) +
  geom_point(
    data = cell_types,
    aes(x = x, y = y, fill = factor(cluster_id)),
    shape = 21, size = 4.2, stroke = 0, colour = "white"
  ) +
  geom_text(
    data = cell_types,
    aes(x = x, y = y, label = cluster_id),
    colour = "white", fontface = "bold", size = 2.05
  ) +
  geom_text(
    data = cell_types,
    aes(x = x + 13, y = y, label = plot_label),
    hjust = 0, parse = TRUE, size = 2.75
  ) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  coord_cartesian(
    xlim = c(0, 780),
    ylim = c(0, 326),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

dir.create("result")
ggsave(filename = "result/nc_forest.pdf", plot = p, width = 7, height = 3, units = "in", device = cairo_pdf, bg = "white")

ggsave(filename = "result/nc_forest.png", plot = p, width = 7, height = 3, units = "in", dpi = 600, bg = "white")
