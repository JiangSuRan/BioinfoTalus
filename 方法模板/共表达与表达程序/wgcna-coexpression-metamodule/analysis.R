## WGCNA 软阈值筛选、动态模块合并与临床性状关联分析
library(WGCNA)

expr_mat <- as.matrix(readRDS("Data/bulk_expr_matrix.rds"))
trait_df <- read.csv("Data/wgcna_sample_traits.csv", row.names = 1, check.names = FALSE)

# 显式对齐表达矩阵样本列与临床性状表行名，防止性状相关性错位
common_samples <- intersect(colnames(expr_mat), rownames(trait_df))
datExpr <- t(expr_mat[, common_samples, drop = FALSE])
trait_df <- trait_df[common_samples, , drop = FALSE]

# 剔除缺失率过高或零方差的基因与离群稀疏样本
gsg <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
trait_df <- trait_df[rownames(datExpr), , drop = FALSE]

# 在 1:20 候选幂次中选取无标度拟合优度 R^2 达 0.85 的最小软阈值
powers <- 1:20
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  RsquaredCut = 0.85,
  verbose = 0
)
softPower <- ifelse(is.na(sft$powerEstimate), 6L, sft$powerEstimate)

# 有向邻接矩阵与拓扑重叠相异度：type = "signed" 区分正负共表达方向
adj_mat <- adjacency(datExpr, type = "signed", power = softPower)
TOM <- TOMsimilarity(adj_mat, TOMType = "signed")
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 基于 TOM 距离的混合动态剪切与近邻模块合并，抑制过度碎片化
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  method = "hybrid",
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)
dynamicColors <- labels2colors(dynamicMods)

merged <- mergeCloseModules(
  exprData = datExpr,
  colors = dynamicColors,
  cutHeight = 0.25,
  verbose = 0
)
moduleColors <- merged$colors
MEs <- orderMEs(merged$newMEs)

# 模块特征基因与临床性状 Pearson 相关及 Student 渐近显著性检验
moduleTraitCor <- cor(MEs, trait_df, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

gene_module <- data.frame(
  gene = colnames(datExpr),
  dynamic_color = dynamicColors,
  merged_color = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(gene_module, "results/wgcna_gene_modules.csv", row.names = FALSE)

wgcna_results <- list(
  softPower = softPower,
  sft = sft,
  geneTree = geneTree,
  dynamicColors = dynamicColors,
  moduleColors = moduleColors,
  MEs = MEs,
  moduleTraitCor = moduleTraitCor,
  moduleTraitPvalue = moduleTraitPvalue,
  gene_module = gene_module
)
saveRDS(wgcna_results, "results/wgcna_modules.rds")
