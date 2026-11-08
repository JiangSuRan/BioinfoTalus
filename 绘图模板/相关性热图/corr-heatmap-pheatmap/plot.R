## 细胞类型比例 Pearson 相关性热图（pheatmap）
library(ggplot2)
library(pheatmap)

rt <- read.csv("Data/cell_type_for_correlation.csv", sep = ",",
               header = TRUE, check.names = FALSE, row.names = 1)
cor_cell_type <- cor(rt)
hm <- pheatmap(cor_cell_type,
               color = colorRampPalette(c("white", "#1874CD", "#FFC125", "#EE4000", "#CD0000"))(1000),
               scale = "none",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               border_color = NA,
               number_color = "black",
               display_numbers = TRUE)
ggsave(filename = "Data/hm_total_correlation.pdf", plot = hm, width = 6.6, height = 6, units = "in")
