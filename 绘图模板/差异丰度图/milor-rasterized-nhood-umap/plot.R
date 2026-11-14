## 大规模图谱 Milo 邻域差异丰度栅格化 UMAP 组合图
library(miloR)
library(scater)
library(ggrastr)
library(ggplot2)
library(patchwork)

milo_obj <- readRDS("Data/milo_object.rds")
da_results <- readRDS("Data/milo_da_results.rds")

milo_obj <- buildNhoodGraph(milo_obj, overlap = 20)

# 点数与边数在十万级以上时栅格化几何层，保留矢量坐标轴与图例避免 PDF 体积膨胀
p_umap <- plotReducedDim(milo_obj, dimred = "UMAP", colour_by = "cell_type", text_by = "cell_type", point_size = 0.5) +
  guides(fill = "none")
p_nhood <- plotNhoodGraphDA(milo_obj, da_results, res_column = "logFC", layout = "UMAP", alpha = 0.1, node_stroke = 0)

p_combined <- rasterize(p_umap, layers = "GeomPoint", dpi = 500) +
  rasterize(p_nhood, layers = c("GeomPoint", "GeomEdgeSegment"), dpi = 500) +
  plot_layout(guides = "collect")

ggsave(filename = "results/milo_rasterized_umap_nhood.pdf", plot = p_combined, width = 20, height = 8, units = "in")
