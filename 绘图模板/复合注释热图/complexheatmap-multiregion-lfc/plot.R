## ComplexHeatmap 跨脑区/跨亚群差异表达 Log2FC 一致性分块热图
library(ComplexHeatmap)
library(circlize)

lfc_mat <- readRDS("data/deg_lfc_matrix.rds")
pval_mat <- readRDS("data/deg_pval_matrix.rds")

# 行按"基因_细胞类型"后缀分块、列按"病理变量:脑区"前缀分块，检验同一基因在不同区域的方向一致性
row_split <- sub(".*_", "", rownames(lfc_mat))
column_split <- sub(":.*", "", colnames(lfc_mat))

ht <- Heatmap(
  lfc_mat,
  name = "logFC",
  col = colorRamp2(c(-0.5, 0, 0.5), c("#2166ac", "white", "#b2182b")),
  row_split = row_split,
  column_split = column_split,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  # 仅在显著格子（p<0.05）标记星号，突出跨区域稳健信号
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (pval_mat[i, j] < 0.05) {
      grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 6))
    }
  }
)

pdf("result/multiregion_deg_lfc_heatmap.pdf", width = 10, height = 8)
draw(ht)
dev.off()
