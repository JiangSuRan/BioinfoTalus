## clustree 多分辨率聚类层级演化树图
library(Seurat)
library(clustree)
library(ggplot2)

seurat_object <- readRDS("Data/choir_multires_seurat.rds")

# CHOIR 层级树（full_tree 含各层 L0/L1/... 聚类列）或 metadata 多分辨率列
ct <- clustree(seurat_object@meta.data, prefix = "L")

ggsave(filename = "results/clustree_resolution_hierarchy.pdf", plot = ct, width = 14, height = 10, units = "in")

# 对照：CHOIR 树对象直接绘制
choir_data <- readRDS("Data/choir_output.rds")
ct_tree <- clustree(choir_data$clusters$full_tree, prefix = "L")
ggsave(filename = "results/clustree_choir_fulltree.pdf", plot = ct_tree, width = 20, height = 10, units = "in")
