## splatter 单细胞 RNA-seq 计数模拟（Splat 方法）
library(splatter)
library(scater)

# 参数估计（可选）：从真实数据估计 Splat 参数后模拟，或直接用默认
params <- newSplatParams(
  nGenes = 1000,
  batchCells = c(500, 500),
  group.prob = c(0.5, 0.3, 0.2)
)

# groups 方法：模拟离散群（细胞类型）；paths 方法模拟连续轨迹
sim <- splatSimulate(params, method = "groups", verbose = FALSE)

# 标准化降维核对模拟结构
sim <- logNormCounts(sim)
sim <- runPCA(sim)
plotPCA(sim, colour_by = "Group")

saveRDS(sim, "results/splatter_sim_groups.rds")
