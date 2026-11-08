## BayesPrism 基于单细胞参考图谱的 Bulk RNA-seq 贝叶斯反卷积
library(BayesPrism)

set.seed(12345)

sc_dat <- readRDS("Data/sc_raw_counts.rds")
sc_meta <- readRDS("Data/sc_metadata.rds")
bk_dat <- readRDS("Data/bulk_counts.rds")

cell.type.labels <- as.character(sc_meta$cell_type)
cell.state.labels <- as.character(sc_meta$cell_state)

# 剔除线粒体/核糖体/重复序列等异常基因群，降低非信息基因对后验估计的干扰
sc_dat.filtered <- cleanup.genes(
  input = sc_dat,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "mito", "rep")
)

# 仅保留蛋白编码基因作为参考特征，与源码引用的 sc.dat.filtered.pc 变量对应
sc_dat.filtered.pc <- select.gene.type(sc_dat.filtered, type = "protein_coding", species = "hs")

myPrism <- new.prism(
  reference = sc_dat.filtered.pc,
  mixture = bk_dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

bp.res <- run.prism(prism = myPrism, n.cores = 8)

# theta 为样本 × 细胞类型的后验比例矩阵（final theta 由 cell state 层级聚合）
theta <- get.fraction(bp = bp.res, which.theta = "final")

saveRDS(bp.res, "results/bayesprism_res.rds")
write.csv(theta, "results/bayesprism_fraction.csv")
