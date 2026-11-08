## ArchR 双细胞评分过滤（addDoubletScores）与迭代 LSI 降维
library(ArchR)

addArchRGenome("hg38")
addArchRThreads(threads = 4)

proj <- readRDS("Data/archr_project.rds")

# TSS 富集与片段数双阈值精筛低质量细胞
proj <- proj[proj$TSSEnrichment > 2 & proj$nFrags > 2500]

# 双细胞先以模拟合并片段的 UMAP 距离评分，再按样本内双细胞率滤除
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

# 迭代 LSI 在 Tile 矩阵上交替聚类-重选特征，抑制测序深度相关的技术维
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),
  varFeatures = 25000,
  dimsToUse = 1:30
)

saveRDS(proj, "results/archr_lsi_project.rds")
