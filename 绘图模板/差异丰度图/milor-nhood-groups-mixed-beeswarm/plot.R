## Milo 混合邻域过滤与邻域社区分组蜂群图
library(miloR)
library(ggplot2)
library(patchwork)

milo_obj <- readRDS("Data/milo_object.rds")
da_results <- readRDS("Data/milo_da_results.rds")

da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "cell_type")
# 主成分比例低于 0.7 的过渡邻域归入 Mixed，防止跨群边缘误判为单一亚型效应
da_results$cell_type <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$cell_type)

p_celltype <- plotDAbeeswarm(da_results, group.by = "cell_type")

set.seed(42)
da_results <- groupNhoods(milo_obj, da_results, max.lfc.delta = 10, overlap = 1)
p_groups <- plotNhoodGroups(milo_obj, da_results, layout = "UMAP") +
  plotDAbeeswarm(da_results, group.by = "NhoodGroup")

ggsave(filename = "results/milo_mixed_filtered_beeswarm.pdf", plot = p_celltype, width = 7, height = 7, units = "in")
ggsave(filename = "results/milo_nhood_groups_panel.pdf", plot = p_groups, width = 14, height = 7, units = "in")
