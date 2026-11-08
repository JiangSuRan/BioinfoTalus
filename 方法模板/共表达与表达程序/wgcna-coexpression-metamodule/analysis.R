## WGCNA 软阈值筛选、拓扑重叠矩阵（TOM）构建与共表达模块识别
library(WGCNA)

datExpr0 <- readRDS("Data/bulk_expr_matrix.rds")
datExpr <- t(datExpr0)

# 无标度拓扑拟合指数达到 0.8 的最小软阈值；无估计时回退默认 6
powers <- c(1:10)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6

# 邻接与拓扑重叠：TOM 缓和单边噪声，signed 网络区分正负相关
adjacency <- adjacency(datExpr, power = softPower, networkType = "signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 动态剪切树：自适应高度切分识别形态各异的模块
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  method = "tree",
  minClusterSize = 30,
  deepSplit = 2
)

# 模块特征基因：每模块表达谱第一主成分，作为模块活性度量
MEs <- moduleEigengenes(datExpr, colors = dynamicMods)$eigengenes
MEs <- orderMEs(MEs)
moduleColors <- labels2colors(dynamicMods)

wgcna_modules <- list(
  softPower = softPower,
  sft = sft,
  geneTree = geneTree,
  modules = dynamicMods,
  moduleColors = moduleColors,
  MEs = MEs
)
saveRDS(wgcna_modules, "results/wgcna_modules.rds")
