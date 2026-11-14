## 样本间亚群比例共现相关性数值热图（ComplexHeatmap）
library(rstatix)
library(ComplexHeatmap)
library(circlize)
library(grid)

freq_df <- readRDS("Data/sample_cluster_freq_matrix.rds")

corr_res <- cor_mat(freq_df)
corr_mat <- as.matrix(corr_res[, -1])
rownames(corr_mat) <- corr_res$rowname

col_fun <- colorRamp2(c(-0.5, 0, 0.5, 1), c("#0F7B9F", "white", "#f4a582", "#D83215"))

ht <- Heatmap(
  corr_mat,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", corr_mat[i, j]), x, y, gp = gpar(fontsize = 6))
  },
  width = ncol(corr_mat) * unit(0.45, "inch"),
  height = nrow(corr_mat) * unit(0.35, "inch"),
  name = "Pearson r"
)

pdf("results/cluster_proportion_cooccurrence_corr.pdf", width = 6, height = 5.5)
draw(ht)
dev.off()
