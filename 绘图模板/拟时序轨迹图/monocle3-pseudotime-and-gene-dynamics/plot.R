## Monocle3 拟时序轨迹图、基因动态拟时曲线与模块热图
library(monocle3)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(stringr)

cds <- readRDS("Data/monocle3_cds.rds")
track_genes <- readRDS("Data/track_genes.rds")

p_ct <- plot_cells(
  cds,
  color_cells_by = "celltype",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 3.5,
  cell_size = 0.7
)

p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = FALSE,
  graph_label_size = 1.5,
  cell_size = 0.7
)

lineage_cds <- cds[rowData(cds)$gene_short_name %in% track_genes, ]
p_dyn <- plot_genes_in_pseudotime(
  lineage_cds,
  color_cells_by = "celltype",
  min_expr = 0.1
)

p_traj <- (p_ct | p_pt) / p_dyn + plot_layout(heights = c(1, 1.1))
ggsave("results/monocle3_trajectory_and_genes.pdf", plot = p_traj, width = 10, height = 8, units = "in")

# 沿主图拓扑检验空间自相关动态基因并按 Leiden 分辨率聚成共表达模块
pr_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
deg_ids <- rownames(subset(pr_res, q_value < 0.05))
mod_df <- find_gene_modules(cds[deg_ids, ])
group_df <- tibble(cell = rownames(colData(cds)), cell_group = as.character(colData(cds)[["celltype"]]))
agg_mat <- aggregate_gene_expression(cds, mod_df, group_df)
rownames(agg_mat) <- str_c("Module ", rownames(agg_mat))

pdf("results/monocle3_gene_modules_heatmap.pdf", width = 6.5, height = 5.5)
pheatmap(agg_mat, cluster_rows = TRUE, cluster_cols = TRUE, scale = "column", clustering_method = "ward.D2", fontsize = 8)
dev.off()
