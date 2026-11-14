## Milo 邻域差异丰度 UMAP 图与亚群蜂群箱线组合图
library(miloR)
library(ggplot2)
library(patchwork)

milo_obj <- readRDS("Data/milo_object.rds")
da_results <- readRDS("Data/milo_da_results.rds")

p_nhood <- plotNhoodGraphDA(milo_obj, da_results, size_range = c(0.8, 3), layout = "UMAP", alpha = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")

da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "cell_type")

p_beeswarm <- plotDAbeeswarm(da_results, group.by = "cell_type", alpha = 1) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  geom_hline(yintercept = 0, linetype = "dashed")

ggsave(filename = "results/milo_nhood_umap.pdf", plot = p_nhood + plot_layout(guides = "collect"), width = 7.5, height = 7.5, units = "in")
ggsave(filename = "results/milo_da_beeswarm.pdf", plot = p_beeswarm, width = 7.5, height = 7.5, units = "in")
