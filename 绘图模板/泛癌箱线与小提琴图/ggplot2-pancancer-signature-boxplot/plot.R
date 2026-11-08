## 泛癌种 T 细胞反应性评分排序箱线散点复合图（ggplot2）
library(ggpubr)
library(dplyr)

res <- read.csv("data/pancancer_patient_score.csv", stringsAsFactors = FALSE)

# 癌种按患者均分升序定级，坐标翻转后自下而上递增排列
mean_score <- res %>%
  group_by(TumorType) %>%
  summarise(mean_score = mean(mean_Tcell_Reactivity1, na.rm = TRUE), .groups = "drop")
res$TumorType <- factor(res$TumorType, levels = mean_score$TumorType[order(mean_score$mean_score)])

tumor_type_colors <- c(
  "BC_adv" = "#9966c7ff", "BC_early" = "#d8b7f9ff", "CC" = "#de5ea3ff",
  "CRC" = "#bd5b60ff", "GBM" = "#e49c6aff", "HCC" = "#8fc887ff",
  "HGSOC" = "#558cc6ff", "HNSCC" = "#f8c352ff", "MEL" = "#a0a1a2ff",
  "NSCLC_adv" = "#389b95ff", "NSCLC_early" = "#84dcd2ff"
)

p <- res %>%
  ggboxplot(x = "TumorType", y = "mean_Tcell_Reactivity1", fill = "TumorType",
            outlier.shape = NA, add = "mean") +
  geom_jitter(size = 0.5, width = 0.3) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  labs(x = NULL, y = NULL, title = "Reactive T-cell Signature") +
  scale_fill_manual(values = tumor_type_colors) +
  coord_flip()

# 仅对 Wilcoxon p<0.05 的癌种对标注显著性
ttests <- compare_means(mean_Tcell_Reactivity1 ~ TumorType, data = res, method = "wilcox") %>%
  filter(p < 0.05)
p <- p + stat_compare_means(
  method = "wilcox",
  label = "p.signif",
  comparisons = lapply(seq_len(nrow(ttests)), function(i) c(ttests$group1[i], ttests$group2[i])),
  tip.length = 0.005,
  vjust = 1.1
)

ggsave(filename = "result/pancancer_reactivity_boxplot.pdf", plot = p, width = 10, height = 5, units = "in")
