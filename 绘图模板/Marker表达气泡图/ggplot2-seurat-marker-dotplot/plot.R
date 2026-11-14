## Marker 基因亚群表达气泡图（DotPlot 统计重绘与极值截断配色）
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(scales)

sc_obj <- readRDS("Data/seurat_obj.rds")
marker_genes <- readRDS("Data/marker_genes.rds") %>% unique()

# 通过 DotPlot 计算缩放平均表达与检出率，再交由 ggplot2 自定义配色边界
base_dot <- DotPlot(sc_obj, features = marker_genes, group.by = "celltype")
dot_df <- base_dot$data

p <- ggplot(dot_df, aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) +
  geom_point() +
  # 压缩超过阈值 3 的极端缩放值，保持中等表达 Marker 的色阶区分度
  scale_color_gradientn(
    colours = viridis(20),
    limits = c(0, 3),
    oob = squish,
    name = "Scaled expr"
  ) +
  scale_size_continuous(range = c(1, 6), name = "Pct. expr (%)") +
  labs(x = NULL, y = NULL) +
  theme_cowplot() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave("results/marker_dotplot.pdf", plot = p, width = 10, height = 5, units = "in")
