## GSVA / ssGSEA 批量转录组基因签名富集评分计算
library(GSVA)
library(data.table)

expr_mat <- as.matrix(readRDS("Data/bulk_tpm_matrix.rds"))

# 基因集列表：两列（签名名、基因符号）展开为命名列表
sigs <- fread("Data/signature_genes.tsv", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(sigs) <- c("cluster", "gene")
inpList <- lapply(unique(sigs$cluster), function(i) sigs$gene[sigs$cluster == i])
names(inpList) <- unique(sigs$cluster)

# TPM 等连续量用 Gaussian 核密度；mx.diff=TRUE 使得分双向有界于 [-1,1]
gsva_es <- gsva(expr_mat, inpList, method = "gsva", kcdf = "Gaussian", mx.diff = TRUE, parallel.sz = 8)

saveRDS(gsva_es, "results/gsva_scores.rds")
