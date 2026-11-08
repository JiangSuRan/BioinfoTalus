## variancePartition (dream) 队列级 Pseudobulk 多组学表型（年龄/性别/个体/细胞类型）方差解释率分解
library(edgeR)
library(variancePartition)
library(BiocParallel)

count_matrix <- read.csv("Data/pseudobulk_counts.csv", row.names = 1, check.names = FALSE)
meta_data <- read.csv("Data/sample_meta.csv", row.names = 1)

# TMM 归一化后以 dream 观测权重进入混合模型
dge <- DGEList(counts = as.matrix(count_matrix))
dge <- calcNormFactors(dge)

param <- SnowParam(8)

# 固定年龄效应；性别/细胞类型/个体为随机截距，分解各方差来源
form <- ~ Age + (1 | sex) + (1 | celltype) + (1 | sample)

vobjDream <- voomWithDreamWeights(dge, form, meta_data, BPPARAM = param)
vp <- fitExtractVarPartModel(vobjDream, form, meta_data, BPPARAM = param)

write.csv(as.data.frame(vp), "results/variance_partition_results.csv")

vp <- sortCols(vp, FUN = median)
p <- plotVarPart(vp)
ggsave(filename = "results/variance_partition_violin.pdf", plot = p, width = 8, height = 5, units = "in")
