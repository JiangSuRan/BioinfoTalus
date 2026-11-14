## GSVA / ssGSEA 批量转录组基因签名富集评分计算
library(GSVA)
library(BiocParallel)
library(data.table)

expr_mat <- as.matrix(readRDS("Data/bulk_tpm_matrix.rds"))

# 基因集列表按签名名分组展开为命名字符向量列表
sigs <- fread("Data/signature_genes.tsv", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(sigs) <- c("cluster", "gene")
inpList <- split(sigs$gene, sigs$cluster)

# GSVA 2.x 参数对象接口；连续表达用 Gaussian 核，单样本排序评分可改用 ssgseaParam()
gsva_param <- gsvaParam(
  exprData = expr_mat,
  geneSets = inpList,
  minSize = 1,
  maxSize = Inf,
  kcdf = "Gaussian",
  maxDiff = TRUE
)
gsva_es <- gsva(gsva_param, verbose = FALSE, BPPARAM = SerialParam(progressbar = FALSE))

saveRDS(gsva_es, "results/gsva_scores.rds")
