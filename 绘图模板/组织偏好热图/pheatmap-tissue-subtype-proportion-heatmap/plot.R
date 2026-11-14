## 多器官亚群比例归一度热图（含谱系注释与分块间隔）
library(pheatmap)
library(grid)

meta <- readRDS("Data/cell_tissue_lineage_metadata.rds")

tab <- table(meta$tissue, meta$subtype)
# 按亚群总细胞数归一化各器官占比，衡量单一亚群在不同器官中的相对分布权重
prop_mat <- sweep(unclass(tab), 2, colSums(tab), "/")

subtype_anno <- unique(meta[, c("subtype", "major")])
rownames(subtype_anno) <- subtype_anno$subtype
subtype_anno <- subtype_anno[, "major", drop = FALSE]

ord <- rownames(subtype_anno)[order(subtype_anno$major)]
prop_mat <- prop_mat[, ord, drop = FALSE]
subtype_anno <- subtype_anno[ord, , drop = FALSE]

gaps <- cumsum(rle(as.character(subtype_anno$major))$lengths)
gaps <- gaps[-length(gaps)]

bk <- c(seq(0, 0.2, by = 0.01), seq(0.21, 0.6, by = 0.01))
cols <- c(
  colorRampPalette(c("white", "firebrick"))(30),
  colorRampPalette(c("firebrick", "Maroon"))(30)
)

p <- pheatmap(
  prop_mat,
  breaks = bk,
  color = cols,
  gaps_col = gaps,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  display_numbers = FALSE,
  annotation_col = subtype_anno,
  border_color = "grey90"
)

pdf("results/tissue_subtype_proportion_heatmap.pdf", width = 14, height = 6)
grid.newpage()
grid.draw(p$gtable)
dev.off()
