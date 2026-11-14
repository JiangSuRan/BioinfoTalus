## 分裂小提琴图展示病例与对照组多基因表达分布及显著性检验
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)

expr_mat <- readRDS("Data/expr_matrix.rds")
sample_meta <- readRDS("Data/sample_group.rds")
target_genes <- readLines("Data/target_genes.txt")

# 将宽表转为长表并按样本对齐分组，便于多基因分面配对比较
plot_df <- expr_mat[intersect(target_genes, rownames(expr_mat)), , drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  melt(id.vars = "Sample", variable.name = "Gene", value.name = "Expression") %>%
  inner_join(sample_meta, by = "Sample") %>%
  mutate(
    Group = factor(Group, levels = c("Control", "Treat")),
    Gene = factor(Gene, levels = intersect(target_genes, rownames(expr_mat)))
  )

# 自定义分裂小提琴几何对象，使对照组与病例组分别占据左右半侧便于对称直观对比
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
    if (length(draw_quantiles) > 0 && !scales::zero_range(range(data$y))) {
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

p <- ggplot(plot_df, aes(x = Gene, y = Expression, fill = Group)) +
  geom_split_violin(trim = FALSE, colour = "grey30", linewidth = 0.4, alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, position = position_dodge(width = 0.22), colour = "grey15", alpha = 0.9) +
  stat_compare_means(aes(group = Group), method = "wilcox.test", label = "p.signif", label.y = max(plot_df$Expression, na.rm = TRUE) * 1.06, hide.ns = FALSE) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Treat" = "#E64B35")) +
  labs(x = "", y = "Normalized Expression", fill = "Group") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "top"
  )

ggsave(filename = "results/split_violin_gene_expression.pdf", plot = p, width = 8, height = 5, units = "in")
