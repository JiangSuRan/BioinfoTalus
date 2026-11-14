## ConsensusClusterPlus 无监督分子分型与 PAC 最优簇数判定
library(ConsensusClusterPlus)

expr_mat <- as.matrix(readRDS("Data/bulk_expr_matrix.rds"))
feature_genes <- readLines("Data/subtyping_genes.txt")
feature_genes <- intersect(feature_genes, rownames(expr_mat))

# 对候选分子标志基因按行中位数中心化，消除基线表达量级差异
expr_sub <- expr_mat[feature_genes, , drop = FALSE]
expr_centered <- sweep(expr_sub, 1, apply(expr_sub, 1, median, na.rm = TRUE))

k_max <- 6L
cc_res <- ConsensusClusterPlus(
  d = expr_centered,
  maxK = k_max,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  distance = "pearson",
  seed = 123,
  plot = "pdf",
  title = "results/consensus_cluster",
  writeTable = FALSE
)

# PAC（模糊聚类比例）量化 CDF 在 [0.1, 0.9] 区间增量，PAC 最小对应最稳定分型 K
k_grid <- 2:k_max
pac_vals <- sapply(k_grid, function(k) {
  cons_mat <- cc_res[[k]]$consensusMatrix
  fn <- ecdf(cons_mat[lower.tri(cons_mat)])
  fn(0.9) - fn(0.1)
})
pac_df <- data.frame(k = k_grid, PAC = pac_vals)
best_k <- pac_df$k[which.min(pac_df$PAC)]

subtype_df <- data.frame(
  Sample = names(cc_res[[best_k]]$consensusClass),
  Subtype = paste0("C", cc_res[[best_k]]$consensusClass),
  BestK = best_k,
  stringsAsFactors = FALSE
)

write.csv(pac_df, "results/consensus_pac_metrics.csv", row.names = FALSE)
write.csv(subtype_df, "results/consensus_subtypes.csv", row.names = FALSE)
saveRDS(cc_res, "results/consensus_results.rds")
