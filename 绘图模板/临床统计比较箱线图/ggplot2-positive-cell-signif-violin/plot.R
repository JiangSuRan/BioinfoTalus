## 阳性表达细胞子集分组小提琴散点与 ggsignif 成对检验图
library(Seurat)
library(ggplot2)
library(ggsignif)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
target_gene <- readRDS("Data/target_gene.rds")[[1]]

expr_vec <- GetAssayData(sc_obj, assay = "RNA", layer = "data")[target_gene, ]
# 剔除高稀疏零表达细胞，避免零膨胀主导小提琴核密度估计及均值检验
pos_df <- sc_obj@meta.data %>%
  mutate(exp = expr_vec) %>%
  filter(exp > 0) %>%
  mutate(group = factor(group))

cmp_list <- combn(levels(pos_df$group), 2, simplify = FALSE)

p <- ggplot(pos_df, aes(x = group, y = exp, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.85) +
  geom_jitter(shape = 16, position = position_jitter(0.18), size = 0.6, alpha = 0.6) +
  geom_signif(
    comparisons = cmp_list,
    step_increase = 0.1,
    map_signif_level = FALSE,
    test = "t.test",
    size = 0.6,
    textsize = 4
  ) +
  labs(title = target_gene, x = "Group", y = "Normalized expression") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, face = "bold", angle = 25, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

ggsave("results/positive_cell_signif_violin.pdf", plot = p, width = 6.5, height = 5.5, units = "in")
