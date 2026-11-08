## SCENIC 单细胞转录因子共表达网络（GENIE3）与 Regulon 活性推断
library(SCENIC)

exprMat <- as.matrix(readRDS("Data/expr_matrix.rds"))
scenicOptions <- initializeScenic(org = "hgnc", dbDir = "Data/scenic", nCores = 10)

genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene = 0.5, minSamples = 3)

# 剔除核糖体/高持家黑名单基因，避免其高共现性产生伪调控子
f.blacklist <- grepl("^RP[LS]", rownames(exprMat), perl = TRUE) | rownames(exprMat) == "MALAT1"
genesKept <- setdiff(genesKept, rownames(exprMat)[f.blacklist])
exprMat_filtered <- exprMat[genesKept, ]

# GENIE3 前先算基因间相关矩阵，仅对高相关对推断树模型以压缩计算量
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)

# 共表达模块 -> motif 验证的 Regulon -> AUCell 逐细胞活性打分
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
aucell <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

saveRDS(aucell, "results/scenic_aucell.rds")
saveRDS(scenicOptions, "results/scenic_options.rds")
