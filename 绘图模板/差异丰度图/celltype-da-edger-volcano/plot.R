## Pseudobulk 亚群差异丰度火山图
library(ggplot2)
library(ggrepel)

da_df <- read.csv("Data/celltype_da_results.csv", stringsAsFactors = FALSE)

p <- ggplot(da_df, aes(x = logFC, y = -log10(FDR), fill = cell_type)) +
  geom_point(pch = 21, size = 5, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linewidth = 0.5, linetype = "dotted", color = "grey50") +
  geom_text_repel(data = da_df[da_df$FDR < 0.1, ], aes(label = cell_type), size = 4, color = "grey30", force = 9, max.overlaps = Inf) +
  labs(x = "Log2(FC)", y = "-Log10(FDR)") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(filename = "results/celltype_da_volcano.pdf", plot = p, width = 8, height = 6, units = "in")
