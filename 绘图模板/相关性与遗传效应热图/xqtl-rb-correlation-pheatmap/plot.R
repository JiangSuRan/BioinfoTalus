## pheatmap 跨细胞类型 cis-eQTL / cis-caQTL 遗传效应相关性（rb）分段断点聚类热图
library(pheatmap)

eQTL_rb_csv <- read.csv("data/eqtl_rb_matrix.csv", row.names = 1, check.names = FALSE)
celltype_sort <- read.csv("data/celltype_sort.csv", stringsAsFactors = FALSE)

# 连字符统一为点号、对角线置 1，并按谱系排序表重排行列
rownames(eQTL_rb_csv) <- gsub("-", ".", rownames(eQTL_rb_csv))
colnames(eQTL_rb_csv) <- gsub("-", ".", colnames(eQTL_rb_csv))
for (ct in rownames(eQTL_rb_csv)) {
  eQTL_rb_csv[ct, ct] <- 1
}
ord <- gsub("-", ".", celltype_sort$celltype[celltype_sort$celltype %in% rownames(eQTL_rb_csv)])
eQTL_rb_csv <- eQTL_rb_csv[ord, ord]

# 谱系注释条：按细胞类型名前缀归并为髓系/B/CD4 T/CD8 T & NK 四系
lineage <- ifelse(grepl("^(ncMono|cMono|intMono|cDC|pDC|AS_DC)", ord), "Myeloid",
           ifelse(grepl("^(Bn|Transitional_B|.*Bm|Plasma|Plasmablast)", ord), "B",
           ifelse(grepl("^CD4", ord), "CD4 T", "CD8 T & NK")))
annotation_col <- data.frame(lineage = factor(lineage), row.names = ord)
ann_colors <- list(lineage = c("Myeloid" = "#E64B35", "B" = "#4DBBD5", "CD4 T" = "#00A087", "CD8 T & NK" = "#3C5488"))

custom_colors <- colorRampPalette(c("#FFC787", "white", "#540D6E"))(100)

# 分段断点：0.8 以下细分中低相关差异，0.8-1.0 段加密放大近对角高相关分辨率
min_val <- min(eQTL_rb_csv)
max_val <- max(eQTL_rb_csv)
breaks_1 <- seq(min_val, 0.8, length.out = 51)
breaks_2 <- seq(0.8, max_val, length.out = 51)
breaks <- c(breaks_1, breaks_2[-1])

pheatmap(
  eQTL_rb_csv,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = custom_colors,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  filename = "results/xqtl_rb_correlation_heatmap.pdf",
  width = 10,
  height = 10
)
