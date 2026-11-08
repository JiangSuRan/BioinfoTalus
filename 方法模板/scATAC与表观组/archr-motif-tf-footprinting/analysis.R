## ArchR 转录因子 Motif 富集与 Tn5 偏倚校正足迹分析（getFootprints）
library(ArchR)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

proj <- readRDS("Data/archr_motif_project.rds")

# 足迹计算依赖按群聚合的 Tn5 插入覆盖率（含正负链偏倚矩阵）
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters2")

# motifPositions 为 motif 名 -> 全基因组结合位点 GRanges；按目标 TF 名模糊匹配展开
motifPositions <- getPositions(proj)
motifs <- c("TBX21", "SPIB", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = proj,
  positions = motifPositions[markerMotifs],
  groupBy = "Clusters2"
)

saveRDS(seFoot, "results/archr_footprints.rds")
