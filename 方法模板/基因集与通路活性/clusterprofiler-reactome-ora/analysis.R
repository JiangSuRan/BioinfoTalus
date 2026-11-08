## clusterProfiler 多基因模块 Reactome / GO 批量过代表富集分析（compareCluster）
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(dplyr)

df <- read.csv("Data/nmf_marker_genes.csv", stringsAsFactors = FALSE)

# 每个模块取前 50 标记基因并映射 ENTREZID 以匹配 ReactomePA 注释
n_top <- 50
list.genes <- lapply(0:11, function(c) {
  x <- df[seq_len(n_top), paste0("gene_NMF", c)] %>% unlist(use.names = FALSE)
  eg <- bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  eg$ENTREZID
})
names(list.genes) <- paste0("NMF", 0:11)

enrich_res <- compareCluster(geneCluster = list.genes, fun = enrichPathway, organism = "human")

write.csv(as.data.frame(enrich_res), "results/reactome_enrichment.csv", row.names = FALSE)
