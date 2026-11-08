## 髓系-T细胞免疫检查点跨临床分组表达热图与点图（DotPlot/pheatmap）
library(pheatmap)
library(RColorBrewer)

mat <- read.csv("data/checkpoint_avg_expr.csv", row.names = 1, check.names = FALSE)
row_anno <- read.csv("data/gene_annotation.csv", row.names = 1, stringsAsFactors = FALSE)
col_anno <- read.csv("data/celltype_annotation.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# 行侧注记基因功能分类（抑制/激活、受体/配体、来源细胞），列侧注记多层级细胞类型与临床亚型
row_anno <- row_anno[rownames(mat), , drop = FALSE]
col_anno <- col_anno[colnames(mat), , drop = FALSE]

mycolors1 <- c("#d73027", "#4575b4")
names(mycolors1) <- unique(row_anno$Mod)
mycolors2 <- c("#542788", "#1b7837")
names(mycolors2) <- unique(row_anno$Flag)
mycolors3 <- c("#7fc97f", "#beaed4", "#c51b7d")
names(mycolors3) <- unique(row_anno$ExpCellType)

mycolors4 <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(col_anno$Major_CellType)))
names(mycolors4) <- unique(col_anno$Major_CellType)
mycolors5 <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(col_anno$Minor_CellType)))
names(mycolors5) <- unique(col_anno$Minor_CellType)
mycolors6 <- c("#b2182b", "#1F618D", "#f4a582", "#9B59B6", "#85929E")
names(mycolors6) <- unique(col_anno$Subset_CellType)
mycolors7 <- c("blue", "pink", "red")
names(mycolors7) <- unique(col_anno$BrCa_type)

anno <- list(
  Mod = mycolors1,
  Flag = mycolors2,
  ExpCellType = mycolors3,
  Major_CellType = mycolors4,
  Minor_CellType = mycolors5,
  Subset_CellType = mycolors6,
  BrCa_type = mycolors7
)

phet <- pheatmap(
  mat,
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  border_color = FALSE,
  fontsize = 5,
  cellheight = 5,
  drop_levels = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 5,
  annotation_row = row_anno,
  annotation_colors = anno,
  annotation_col = col_anno,
  cellwidth = 5
)

pdf("result/checkpoint_clinical_heatmap.pdf", width = 10, height = 10)
print(phet)
dev.off()
