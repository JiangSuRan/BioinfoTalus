## 肿瘤微环境亚型与多维临床特征整合景观热图（ComplexHeatmap）
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

intercells <- read.csv("data/microenvironment_mat.csv", row.names = 1, check.names = FALSE)
subtype_anno <- read.csv("data/subtype_annotation.csv", row.names = 1, check.names = FALSE)

subtype_anno <- subtype_anno[colnames(intercells), , drop = FALSE]
subtype_cols <- c("NT" = "#8c96c6", "PT" = "#8856a7", "HM" = "#9e9ac8")
ann_colors <- list(subtype = subtype_cols[unique(subtype_anno$subtype)])

# 互作计数 log1p 后按列 scale（不中心化）并截断到 2，突出高互作组合同时保留动态范围
pltd <- apply(log1p(intercells), 2, scale, center = FALSE)
rownames(pltd) <- rownames(intercells)
pltd[pltd > 2] <- 2

p <- ComplexHeatmap::pheatmap(
  pltd,
  cluster_cols = FALSE,
  annotation_col = subtype_anno,
  annotation_colors = ann_colors,
  cellwidth = 10,
  cellheight = 10,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)
)

pdf("result/pdac_microenvironment_landscape.pdf", width = 10, height = 8)
draw(p)
dev.off()
