## 关键基因与免疫浸润细胞 Spearman 相关性方形气泡热图
library(dplyr)
library(reshape2)
library(ggplot2)

gene_expr <- readRDS("Data/hub_gene_expression.rds")
immune_mat <- readRDS("Data/immune_infiltration_matrix.rds")

common_samples <- intersect(colnames(gene_expr), colnames(immune_mat))
gene_sub <- t(gene_expr[, common_samples, drop = FALSE])
immune_sub <- t(immune_mat[, common_samples, drop = FALSE])

# 逐对计算 Spearman 相关系数及显著性 P 值，适应非正态免疫评分分布
cor_list <- vector("list", ncol(gene_sub) * ncol(immune_sub))
idx <- 1L
for (g in colnames(gene_sub)) {
  for (im in colnames(immune_sub)) {
    res <- suppressWarnings(cor.test(as.numeric(gene_sub[, g]), as.numeric(immune_sub[, im]), method = "spearman", exact = FALSE))
    cor_list[[idx]] <- data.frame(
      Gene = g,
      ImmuneCell = im,
      rho = unname(res$estimate),
      pvalue = res$p.value,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

cor_df <- bind_rows(cor_list) %>%
  mutate(
    abs_rho = abs(rho),
    sig_label = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01 ~ "**",
      pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Gene = factor(Gene, levels = rev(colnames(gene_sub))),
    ImmuneCell = factor(ImmuneCell, levels = colnames(immune_sub))
  )

# 采用空心可填充方形点（shape = 22），同时编码相关方向与强度并叠加星号显著性
p <- ggplot(cor_df, aes(x = ImmuneCell, y = Gene)) +
  geom_point(aes(size = abs_rho, fill = rho), shape = 22, colour = "grey25", stroke = 0.5) +
  geom_text(aes(label = sig_label), colour = "black", size = 3.5, vjust = 0.72) +
  scale_fill_gradient2(low = "#2B5C8F", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  scale_size_continuous(range = c(2.5, 9.5), limits = c(0, 1), name = "|rho|") +
  labs(x = "", y = "") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black")
  )

ggsave(filename = "results/gene_immune_square_bubble_heatmap.pdf", plot = p, width = 10, height = 4.8, units = "in")
