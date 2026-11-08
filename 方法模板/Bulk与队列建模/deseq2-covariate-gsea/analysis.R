## DESeq2 多协变量差异表达拟合与排序基因列表 GSEA 通路分析
library(DESeq2)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

df <- read.delim("Data/bulk_counts.tsv", stringsAsFactors = FALSE, check.names = FALSE)

# 同名基因计数求和并以基因符号建行名
df <- df %>%
  select(-Gene_id) %>%
  group_by(Gene_name) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("Gene_name")

df.cli <- read.delim("Data/clinical_metadata.tsv", stringsAsFactors = FALSE) %>%
  filter(id %in% colnames(df))

dds <- DESeqDataSetFromMatrix(countData = df[, df.cli$id], colData = df.cli, design = ~ disease)
dds$disease <- relevel(dds$disease, ref = "HC")

# 低计数过滤：阈值计数须在至少最小组样本数中达到
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds)

res.tidy <- results(dds, tidy = TRUE)
write.csv(res.tidy, "results/deseq2_results.csv", row.names = FALSE)

# 按 log2FC 排序的基因列表做 Reactome GSEA
gL <- res.tidy[, c("row", "log2FoldChange")] %>%
  arrange(desc(log2FoldChange)) %>%
  as.data.frame()
eg <- bitr(gL$row, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gL <- gL %>%
  left_join(eg, by = c("row" = "SYMBOL")) %>%
  drop_na() %>%
  column_to_rownames("ENTREZID")
geneList <- gL$log2FoldChange
names(geneList) <- rownames(gL)

gsea_res <- gsePathway(geneList, organism = "human")
write.csv(as.data.frame(gsea_res), "results/gsea_results.csv", row.names = FALSE)
