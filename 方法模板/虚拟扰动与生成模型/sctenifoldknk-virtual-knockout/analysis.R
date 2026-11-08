## scTenifoldKnk 单细胞基因调控网络构建与单基因虚拟敲除（Virtual KO）差异调控检验
library(scTenifoldKnk)

count_mat <- readRDS("Data/count_matrix.rds")

# 500 细胞子采样 x 10 次主成分回归网络 + 张量分解去噪；KO 网络经流形对齐与 WT 比较并做差异调控检验
res <- scTenifoldKnk(
  countMatrix = count_mat,
  gKO = "TargetGene",
  qc = TRUE,
  nc_nNet = 10,
  nc_nCells = 500
)

saveRDS(res, "results/sctenifoldknk_res.rds")
write.csv(res$diffRegulation, "results/virtual_ko_diff_regulation.csv", row.names = FALSE)
