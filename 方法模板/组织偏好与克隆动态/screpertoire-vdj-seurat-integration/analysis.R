## scRepertoire 多样本 VDJ 配对整合、克隆回写与 STARTRAC / 克隆偏好定量
library(scRepertoire)
library(Seurat)
library(dplyr)

contig_list <- readRDS("Data/contig_list.rds")
sc_obj <- readRDS("Data/seurat_obj.rds")

combined_tcr <- combineTCR(contig_list, samples = names(contig_list))

sc_vdj <- combineExpression(
  combined_tcr,
  sc_obj,
  cloneCall = "gene",
  group.by = "sample",
  proportion = FALSE,
  cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500)
)

# 量化各亚群克隆扩增、跨组织迁移与状态转换 STARTRAC 指数及亚群偏倚
startrac_tbl <- StartracDiversity(sc_vdj, cloneCall = "aa", type = "orig.ident", group.by = "seurat_clusters", exportTable = TRUE)
bias_tbl <- clonalBias(sc_vdj, cloneCall = "aa", split.by = "orig.ident", group.by = "seurat_clusters", n.boots = 10, min.expand = 5, exportTable = TRUE)

saveRDS(combined_tcr, "results/combined_tcr.rds")
saveRDS(sc_vdj, "results/seurat_vdj.rds")
write.table(startrac_tbl, "results/startrac_diversity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bias_tbl, "results/clonal_bias.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
