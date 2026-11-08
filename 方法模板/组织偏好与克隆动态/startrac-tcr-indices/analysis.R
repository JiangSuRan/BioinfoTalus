## STARTRAC 单细胞 TCR 克隆扩增、跨组织迁移与状态转换指数计算
library(Startrac)

in_dat <- readRDS("Data/tcr_cell_metadata.rds")

# STARTRAC 以 majorCluster 与 clone.id 为亚群和克隆标识列名，先统一字段
in_dat$majorCluster <- as.character(in_dat$meta.cluster)
in_dat$clone.id <- in_dat$cloneID

# patient × majorCluster 细胞数不足 10 的稀疏组合无法支撑克隆结构估计，先剔除
ncell.patient.cluster <- sort(unclass(in_dat[, table(sprintf("%s.%s", patient, majorCluster))]))
in_dat <- in_dat[ncell.patient.cluster[sprintf("%s.%s", patient, majorCluster)] >= 10, ]

# n.perm 置换次数决定指数显著性检验的最小可分辨 p 值（1/(n.perm+1)）
out <- Startrac.run(in_dat, proj = "panC", verbose = FALSE, cores = 4, n.perm = 1000)

saveRDS(out, "results/startrac_output.rds")
