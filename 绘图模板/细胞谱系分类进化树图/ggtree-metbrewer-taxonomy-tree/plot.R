## ggtree 聚类树 + 堆叠小提琴组合图（树图对齐细胞簇标签）
library(treeio)
library(ggtree)
library(ape)
library(Seurat)
library(patchwork)

seu <- readRDS("Data/subset_seurat.rds")
Idents(seu) <- "region_cluster"

# 聚类均值表达谱 → 层次聚类 → ggtree 树图（tiplab 对齐簇名）
cl_means <- sapply(split(seu, Idents(seu)), function(x) {
  rowMeans(GetAssayData(x, assay = "RNA", layer = "counts")[VariableFeatures(x), ])
})
hc <- hclust(dist(t(cl_means)), method = "average")
dend <- as.phylo(hc)

p_tree <- ggtree(dend) + scale_x_reverse() + geom_tiplab(size = 3)

# 堆叠小提琴（VlnPlot stack）与树图拼接
p_vln <- VlnPlot(seu, features = VariableFeatures(seu)[1:10], stack = TRUE, flip = TRUE) +
  theme(legend.position = "none")

p <- p_vln + p_tree + patchwork::plot_layout(widths = c(0.7, 0.3))

ggsave(filename = "results/ggtree_vln_combo.pdf", plot = p, width = 10, height = 8, units = "in")
