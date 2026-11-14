## Seurat v5 Xenium 亚细胞空间细胞类型、转录本分子定位与特征表达联合图
library(Seurat)
library(ggplot2)
library(patchwork)

seu <- readRDS("Data/xenium_seurat.rds")

# shape = 21 点边界在 R PDF 设备中传 NA 会整层隐身，须显式指定细描边
p_dim <- ImageDimPlot(
  object = seu,
  fov = "fov",
  group.by = "cell_type",
  molecules = c("EPCAM", "CD3D", "COL1A1"),
  nmols = 5000,
  size = 1.8,
  alpha = 0.75,
  mols.size = 0.35,
  mols.cols = c("#FFBE0B", "#00F5D4", "#FF006E"),
  dark.background = TRUE,
  border.color = "grey15",
  border.size = 0.15
) +
  ggtitle("Cell Identities & Subcellular Transcripts") +
  theme(
    plot.title = element_text(size = 11, face = "bold", color = "white"),
    plot.background = element_rect(fill = "black", color = NA)
  )

# 同视场空间连续基因表达强度分布，q95 截断减弱极高表达像元拉偏色带
p_feat <- ImageFeaturePlot(
  object = seu,
  fov = "fov",
  features = "EPCAM",
  size = 1.8,
  max.cutoff = "q95",
  cols = c("#1B263B", "#F4A261", "#E63946"),
  dark.background = TRUE,
  border.color = "grey15",
  border.size = 0.15
) +
  ggtitle("EPCAM Spatial Expression") +
  theme(
    plot.title = element_text(size = 11, face = "bold", color = "white"),
    plot.background = element_rect(fill = "black", color = NA)
  )

p_combined <- (p_dim + p_feat + plot_layout(ncol = 2)) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "black", color = NA)))
ggsave("results/xenium_spatial_panels.pdf", p_combined, width = 12, height = 5.5, units = "in")
