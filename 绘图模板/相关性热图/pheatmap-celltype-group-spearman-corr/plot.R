## 亚群×实验分组平均表达谱 Spearman 相关性聚类热图
library(pheatmap)
library(viridis)
library(grid)

avg_expr <- readRDS("Data/celltype_group_avg_expr.rds")

cor_mat <- cor(avg_expr, method = "spearman", use = "pairwise.complete.obs")

p <- pheatmap(
  cor_mat,
  color = viridis(256, option = "D"),
  border_color = NA
)

pdf("results/celltype_group_spearman_corr_heatmap.pdf", width = 7, height = 6.5)
grid.newpage()
grid.draw(p$gtable)
dev.off()
