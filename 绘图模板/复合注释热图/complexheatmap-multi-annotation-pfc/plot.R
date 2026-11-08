## ComplexHeatmap 多区域临床表型与细胞状态复合注释热图
library(ComplexHeatmap)
library(circlize)

mat <- readRDS("data/ast_reactivity_mat.rds")
meta <- read.csv("data/sample_annotation.csv", row.names = 1)

# 供体级病理（nrad）与认知（cogdxad）组合为 delta 四分类，与矩阵列（供体）按行名对齐
meta <- meta[colnames(mat), , drop = FALSE]
meta$nrad <- ifelse(meta$niareagansc > 2, "CTRL", "AD")
meta$cogdxad <- ifelse(meta$cogdx < 4, "CTRL", "AD")
meta$delta <- paste(ifelse(meta$nrad == "AD", "AD", "N"), ifelse(meta$cogdxad == "AD", "CI", "N"), sep = "-")

# 供体方向 z-score，使不同基因/评分在同一色标下可比
mat_scaled <- t(scale(t(mat)))

row_anno <- rowAnnotation(
  Signature = ifelse(rownames(mat) %in% c("GFAP", "OSMR", "CD44"), "Reactive", "Glu-uptake"),
  col = list(Signature = c("Reactive" = "#b2182b", "Glu-uptake" = "#2166ac"))
)
top_anno <- HeatmapAnnotation(
  delta = meta$delta,
  nrad = meta$nrad,
  cogdxad = meta$cogdxad,
  col = list(
    delta = c("AD-CI" = "#b2182b", "AD-N" = "#ef8a62", "N-CI" = "#67a9cf", "N-N" = "#2166ac"),
    nrad = c("AD" = "#d73027", "CTRL" = "#4575b4"),
    cogdxad = c("AD" = "#d73027", "CTRL" = "#4575b4")
  )
)

ht <- Heatmap(
  mat_scaled,
  name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
  top_annotation = top_anno,
  left_annotation = row_anno,
  cluster_columns = FALSE,
  show_column_names = FALSE
)

pdf("result/ast_reactivity_heatmap.pdf", width = 8, height = 6)
draw(ht)
dev.off()
