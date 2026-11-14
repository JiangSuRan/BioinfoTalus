## 单细胞水平亚群排序 Marker 表达热图（截断 Z-score 与动态行列间隙）
library(Seurat)
library(dplyr)
library(pheatmap)
library(viridis)
library(scales)

sc_obj <- readRDS("Data/seurat_obj.rds")
marker_df <- readRDS("Data/top_markers_df.rds") %>%
  arrange(cluster) %>%
  distinct(gene, .keep_all = TRUE)

Idents(sc_obj) <- "celltype"
sc_sub <- subset(sc_obj, downsample = 100)
sc_sub <- ScaleData(sc_sub, features = marker_df$gene, verbose = FALSE)

col_anno <- sc_sub@meta.data[, "celltype", drop = FALSE] %>%
  arrange(celltype)

expr_mat <- GetAssayData(sc_sub, assay = "RNA", layer = "scale.data")[marker_df$gene, rownames(col_anno), drop = FALSE]
# 对 Z-score 设定不对称上下界 [-1, 2.5]，避免极个别细胞强信号使热图整体过暗
expr_mat[expr_mat >= 2.5] <- 2.5
expr_mat[expr_mat < -1.5] <- -1

celltypes <- levels(factor(col_anno$celltype))
type_cols <- setNames(hue_pal()(length(celltypes)), celltypes)
gaps_r <- marker_df$cluster %>% factor(levels = celltypes) %>% table() %>% cumsum() %>% head(-1) %>% as.numeric()
gaps_c <- col_anno$celltype %>% factor(levels = celltypes) %>% table() %>% cumsum() %>% head(-1) %>% as.numeric()

pdf("results/cell_marker_grouped_heatmap.pdf", width = 9, height = 7)
pheatmap(
  expr_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = col_anno,
  annotation_colors = list(celltype = type_cols),
  gaps_row = gaps_r,
  gaps_col = gaps_c,
  color = viridis(256, option = "D"),
  fontsize_row = 7
)
dev.off()
