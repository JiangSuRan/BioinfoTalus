## ggstatsplot 自动内嵌统计检验的分组组间比较图
library(ggstatsplot)
library(dplyr)

ggplotdata <- read.csv("Data/sample_fractions.csv", stringsAsFactors = FALSE)
ggplotdata$TissueType <- factor(ggplotdata$TissueType, levels = c("H", "U", "P", "M"))

# grouped_ggbetweenstats：分组小提琴+均值点，自动内嵌组间检验与效应量，FDR 校正
p <- ggplotdata %>%
  ggstatsplot::grouped_ggbetweenstats(
    data = .,
    x = TissueType,
    y = frac,
    grouping.var = g_seurat_clusters,
    xlab = "Tissue type",
    ylab = "Sample fraction",
    pairwise.display = "significant",
    p.adjust.method = "fdr",
    ggtheme = ggstatsplot::theme_ggstatsplot(),
    package = "ggsci",
    palette = "default_jco",
    plotgrid.args = list(ncol = 3)
  )

ggsave(filename = "results/ggstatsplot_grouped_comparison.pdf", plot = p, width = 12, height = 10, units = "in")
