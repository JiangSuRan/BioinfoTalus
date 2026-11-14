## Slingshot 降维多谱系曲线拟合与拟时序亚群密度分布联图
library(SingleCellExperiment)
library(slingshot)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

sce <- readRDS("Data/slingshot_sce.rds")
cluster_col <- "celltype"
rd_name <- "UMAP"

coords <- as.data.frame(reducedDim(sce, rd_name)[, 1:2])
colnames(coords) <- c("Dim1", "Dim2")
coords$Cluster <- factor(colData(sce)[[cluster_col]])

# 提取 Slingshot 平滑主曲线坐标，叠加于单细胞降维图上展示分支分化路径
curves <- slingCurves(sce, as.df = TRUE)
colnames(curves)[1:2] <- c("Dim1", "Dim2")
endpoints <- curves %>%
  group_by(Lineage) %>%
  slice_tail(n = 1) %>%
  ungroup()

n_clus <- nlevels(coords$Cluster)
pal <- colorRampPalette(brewer.pal(min(9, max(3, n_clus)), "Set1"))(n_clus)

p_curve <- ggplot(coords, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(colour = Cluster), size = 0.7, alpha = 0.75) +
  geom_path(data = curves, aes(x = Dim1, y = Dim2, group = Lineage), colour = "black", linewidth = 0.9) +
  geom_label(data = endpoints, aes(x = Dim1, y = Dim2, label = paste0("L", Lineage)), size = 3.2, fontface = "bold", label.size = 0.25) +
  scale_colour_manual(values = pal) +
  labs(title = "Slingshot Lineages", x = paste0(rd_name, "_1"), y = paste0(rd_name, "_2")) +
  theme_classic(base_size = 11) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

pt_df <- as.data.frame(slingPseudotime(sce)) %>%
  mutate(Cluster = factor(colData(sce)[[cluster_col]])) %>%
  melt(id.vars = "Cluster", variable.name = "Lineage", value.name = "Pseudotime") %>%
  filter(!is.na(Pseudotime))

p_density <- ggplot(pt_df, aes(x = Pseudotime, fill = Cluster, colour = Cluster)) +
  geom_density(alpha = 0.35, linewidth = 0.5) +
  facet_wrap(~Lineage, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(title = "Pseudotime Density Along Lineages", x = "Pseudotime", y = "Density") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "none")

p_all <- p_curve + p_density + plot_layout(widths = c(1.25, 1))

ggsave(filename = "results/slingshot_curves_and_density.pdf", plot = p_all, width = 12, height = 5.5, units = "in")
