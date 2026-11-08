## Seurat 空间转录组分区特征打分与空间域差异表达分析
library(Seurat)
library(dplyr)

st_seu <- readRDS("Data/spatial_seurat.rds")

# 线粒体比例上限取 median+2MAD 与 Q3+1.5IQR 的较小者，且以 30% 兜底，避免阈值过紧或过松
st_seu$percent.mt <- PercentageFeatureSet(st_seu, pattern = "^MT-")
med_cutoff <- median(st_seu$percent.mt) + 2 * mad(st_seu$percent.mt)
IQR_cutoff <- as.numeric(quantile(st_seu$percent.mt)[4] + IQR(st_seu$percent.mt) * 1.5)
mitofilter <- max(min(med_cutoff, IQR_cutoff), 30)

st_seu <- subset(st_seu, subset = nCount_Spatial >= 500 & nFeature_Spatial > 500 & percent.mt < mitofilter)

# 空间域标记基因：按分区标识列做 one-vs-rest Wilcoxon 差异检验
Idents(st_seu) <- "region"
region_markers <- FindAllMarkers(st_seu, test.use = "wilcox", only.pos = TRUE)

region_markers <- region_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()

write.csv(region_markers, "results/spatial_region_markers.csv", row.names = FALSE)
