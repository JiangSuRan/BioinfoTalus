## 三层同心圆 ceRNA（lncRNA-miRNA-mRNA）调控网络 ggraph 图
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)

nodes_df <- readRDS("Data/cerna_nodes.rds")
edges_df <- readRDS("Data/cerna_edges.rds")

g <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)

# 按分子层级分配固定同心圆半径并分环反向对齐标注，避免中间环 miRNA 标签向外侵入外环 lncRNA
ring_radii <- c("mRNA" = 1.1, "miRNA" = 3.1, "lncRNA" = 4.2)
coords <- matrix(0, nrow = nrow(nodes_df), ncol = 2)
for (type_name in names(ring_radii)) {
  idx <- which(nodes_df$Type == type_name)
  theta <- seq(0, 2 * pi, length.out = length(idx) + 1L)[seq_along(idx)]
  r <- ring_radii[[type_name]]
  coords[idx, 1] <- r * cos(theta)
  coords[idx, 2] <- r * sin(theta)
}

layout_df <- create_layout(g, layout = "manual", x = coords[, 1], y = coords[, 2]) %>%
  mutate(
    deg_rad = atan2(y, x),
    text_angle = ifelse(Type == "mRNA", 0, ifelse(cos(deg_rad) < 0, deg_rad * 180 / pi + 180, deg_rad * 180 / pi)),
    text_hjust = case_when(
      Type == "mRNA" ~ 0.5,
      Type == "miRNA" ~ ifelse(cos(deg_rad) < 0, -0.18, 1.18),
      TRUE ~ ifelse(cos(deg_rad) < 0, 1.18, -0.18)
    )
  )

p <- ggraph(layout_df) +
  geom_edge_link(colour = "grey72", alpha = 0.55, linewidth = 0.45) +
  geom_node_point(aes(fill = Type, shape = Type, size = Type), colour = "grey20", stroke = 0.5) +
  geom_node_text(
    aes(label = name, angle = text_angle, hjust = text_hjust),
    data = function(df) filter(df, Type != "mRNA"),
    size = 2.8,
    colour = "black"
  ) +
  geom_node_label(
    aes(label = name),
    data = function(df) filter(df, Type == "mRNA"),
    size = 3.0,
    fontface = "bold",
    fill = alpha("white", 0.85),
    label.size = 0,
    label.padding = unit(0.14, "lines")
  ) +
  scale_fill_manual(values = c("mRNA" = "#E64B35", "miRNA" = "#4DBBD5", "lncRNA" = "#00A087")) +
  scale_shape_manual(values = c("mRNA" = 21, "miRNA" = 22, "lncRNA" = 24)) +
  scale_size_manual(values = c("mRNA" = 6.0, "miRNA" = 4.2, "lncRNA" = 3.5)) +
  coord_fixed(xlim = c(-5.8, 5.8), ylim = c(-5.8, 5.8)) +
  labs(fill = "Molecule Type", shape = "Molecule Type", size = "Molecule Type") +
  theme_void(base_size = 11) +
  theme(legend.position = "right")

ggsave(filename = "results/cerna_concentric_ring_network.pdf", plot = p, width = 9.0, height = 8.0, units = "in")
