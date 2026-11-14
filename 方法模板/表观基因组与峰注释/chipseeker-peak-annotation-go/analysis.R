## ChIPseeker 多样本染色质峰基因组注释与 TSS 近端基因 GO 富集
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)

samples_df <- read.table("Data/peak_samples.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
peak_files <- setNames(samples_df$BedPath, samples_df$Sample)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 启动子窗口设定为 TSS 上下游 3 kb，统一注释基因符号与基因组区域分布
peak_anno_list <- lapply(
  peak_files,
  annotatePeak,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db",
  verbose = FALSE
)

anno_tables <- lapply(names(peak_anno_list), function(s_name) {
  df <- as.data.frame(peak_anno_list[[s_name]])
  df$Sample <- s_name
  df
})
merged_anno <- bind_rows(anno_tables)

# 仅保留距 TSS 3 kb 内峰关联基因作为调控候选集，并按样本比较 GO 生物过程
tss_gene_list <- lapply(anno_tables, function(df) {
  df %>%
    filter(!is.na(geneId), abs(distanceToTSS) <= 3000) %>%
    pull(geneId) %>%
    unique()
})
names(tss_gene_list) <- names(peak_anno_list)

go_comp <- compareCluster(
  geneClusters = tss_gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

saveRDS(peak_anno_list, "results/chipseeker_peak_annotations.rds")
write.csv(merged_anno, "results/chipseeker_peak_annotations.csv", row.names = FALSE)
write.csv(as.data.frame(go_comp), "results/chipseeker_tss_go_enrichment.csv", row.names = FALSE)
