## ChIPseeker 多样本染色质峰基因组区段分布、TSS 距离谱与交集 UpSet 图
library(ChIPseeker)
library(ggupset)
library(ggplot2)
library(patchwork)

peak_anno_list <- readRDS("results/chipseeker_peak_annotations.rds")

# 多样本峰注释区段构成对比与相对 TSS 位置富集曲线组合图
p_bar <- plotAnnoBar(peak_anno_list) +
  theme_bw(base_size = 11) +
  ggtitle("Genomic Feature Distribution")

p_tss <- plotDistToTSS(peak_anno_list) +
  theme_bw(base_size = 11) +
  ggtitle("Distance to Transcription Start Site")

p_panels <- (p_bar / p_tss) + plot_layout(heights = c(1, 1))
ggsave("results/chipseeker_multi_sample_panels.pdf", p_panels, width = 9.5, height = 7.5, units = "in")

# 代表样本基因特征重叠 UpSet 组合与 Venn 饼图
pdf("results/chipseeker_upset_vennpie.pdf", width = 8.5, height = 5.5)
print(upsetplot(peak_anno_list[[1]], vennpie = TRUE))
dev.off()
