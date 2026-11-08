## 跨癌种与跨细胞类型标记基因平均表达分面热图
library(ggplot2)
library(dplyr)

expr.df <- read.csv("data/pancancer_avg_expr.csv", stringsAsFactors = FALSE)

genes <- c("PDCD1", "CD274", "PDCD1LG2")
order.tumortype <- rev(c("BC_early", "BC_adv", "CC", "CRC", "GBM", "HCC", "HGSOC", "HNSCC", "MEL", "NSCLC_early", "NSCLC_adv"))
order.celltype <- c("T cell", "NK", "B cell", "Macro/Mono", "DC", "Mast cell", "Cancer/epit cell", "Fibroblast", "EC")

# 基因 × 癌种内跨细胞类型 z-score（分母加 0.01 防零方差），行内对比该癌种各细胞类型相对高低
expr.df <- expr.df %>%
  group_by(Gene, tumortype) %>%
  mutate(zscore = (score - mean(score)) / (sd(score) + 0.01)) %>%
  ungroup()

expr.df$Gene <- factor(expr.df$Gene, levels = genes)
expr.df$celltype <- factor(expr.df$celltype, levels = order.celltype)
expr.df$tumortype <- factor(expr.df$tumortype, levels = order.tumortype)

p <- ggplot(data = expr.df) +
  geom_tile(aes(x = celltype, y = tumortype, fill = zscore, width = 0.92, height = 0.92)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  facet_grid(rows = vars(Gene), scales = "free", space = "free") +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 12),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = NA),
    panel.background = element_rect(fill = NA, color = "black", linetype = "solid"),
    panel.border = element_rect(fill = NA)
  )

ggsave(filename = "result/pancancer_celltype_heatmap.pdf", plot = p, width = 5, height = 8.5, units = "in")
