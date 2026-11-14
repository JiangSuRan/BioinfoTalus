## Marker 基因紧凑堆叠小提琴图（单基因上限刻度与零冗余排布）
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
marker_genes <- readRDS("Data/marker_genes.rds") %>% unique()
Idents(sc_obj) <- "celltype"

# 单基因精简面板：隐藏冗余坐标轴标题，仅保留左置水平基因名以压缩纵向高度
make_vln_panel <- function(obj, gene) {
  p <- VlnPlot(obj, features = gene, pt.size = 0) +
    labs(x = NULL, y = gene, title = NULL) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 9, angle = 0, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 8),
      plot.margin = unit(c(-0.15, 0.2, -0.15, 0.2), "cm")
    )
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range) %>% ceiling()
  p + scale_y_continuous(breaks = ymax) + expand_limits(y = ymax)
}

plot_list <- lapply(marker_genes, make_vln_panel, obj = sc_obj)
last_idx <- length(plot_list)
plot_list[[last_idx]] <- plot_list[[last_idx]] +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.ticks.x = element_line(),
    plot.margin = unit(c(-0.15, 0.2, 0.2, 0.2), "cm")
  )

p_stack <- wrap_plots(plotlist = plot_list, ncol = 1)
ggsave("results/marker_stacked_violin.pdf", plot = p_stack, width = 7, height = max(4, 0.65 * last_idx + 1.2), units = "in")
