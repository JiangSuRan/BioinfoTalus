## SCENIC 特异性 Regulon RSS 排名散点与分组活性热图
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)

rss_mat <- readRDS("Data/rss_matrix.rds") %>% as.data.frame()
auc_mat <- readRDS("Data/regulon_auc.rds") %>% as.matrix()
cell_meta <- readRDS("Data/cell_meta.rds")

make_rss_plot <- function(ct) {
  rank_df <- data.frame(TF = rownames(rss_mat), score = rss_mat[[ct]]) %>%
    filter(!is.na(score)) %>%
    arrange(desc(score)) %>%
    mutate(rank = seq_len(n()))
  top_df <- rank_df %>% slice_head(n = min(6, nrow(rank_df)))

  ggplot(rank_df, aes(x = rank, y = score)) +
    geom_point(size = 1.2, shape = 16, color = "#1F77B4", alpha = 0.45) +
    geom_point(data = top_df, size = 1.8, color = "#DC050C") +
    geom_text_repel(
      data = top_df,
      aes(label = TF),
      color = "black",
      size = 3,
      fontface = "italic",
      arrow = arrow(ends = "first", length = unit(0.01, "npc")),
      box.padding = 0.25,
      point.padding = 0.3,
      segment.size = 0.3,
      force = 1,
      max.iter = 3000
    ) +
    labs(x = "Regulons Rank", y = "Specificity Score", title = ct) +
    theme_bw() +
    theme(
      axis.title = element_text(colour = "black", size = 11),
      axis.text.y = element_text(colour = "black", size = 9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

rss_plots <- lapply(colnames(rss_mat), make_rss_plot)
p_rss <- plot_grid(plotlist = rss_plots, ncol = min(3, length(rss_plots)))
ggsave("results/scenic_rss_rank_panels.pdf", plot = p_rss, width = 10, height = 3.8 * ceiling(length(rss_plots) / 3), units = "in")

shared_cells <- intersect(colnames(auc_mat), rownames(cell_meta))
auc_mat <- auc_mat[, shared_cells, drop = FALSE]
cell_groups <- split(shared_cells, cell_meta[shared_cells, "celltype"])
mean_auc <- sapply(cell_groups, function(cells) rowMeans(auc_mat[, cells, drop = FALSE]))
# 剔除跨亚群无变化的常数值 Regulon，防止行标准差为 0 产生 NaN
keep_reg <- apply(mean_auc, 1, sd) > 0
scaled_auc <- t(scale(t(mean_auc[keep_reg, , drop = FALSE]), center = TRUE, scale = TRUE))
top_n <- min(80, nrow(scaled_auc))

pdf("results/scenic_regulon_activity_heatmap.pdf", width = 7, height = 8)
draw(Heatmap(scaled_auc[seq_len(top_n), , drop = FALSE], name = "Regulon activity", row_names_gp = gpar(fontsize = 7)))
dev.off()
