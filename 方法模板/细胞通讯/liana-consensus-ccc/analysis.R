## LIANA+ 单细胞多算法共识配体-受体通讯推断与聚合
library(liana)

seu <- readRDS("Data/seurat_annotated.rds")

# liana_wrap 一次运行 NATMI/CellPhoneDB/SingleCellSignalNet 等多个方法，以同一配体-受体资源为基准
Idents(seu) <- "cell_type"
liana_res <- liana_wrap(seu, idents_col = "cell_type")

# magnitude 聚合保留连续通讯强度指标，跨方法标准化后给出共识排名
liana_res <- liana_aggregate(liana_res, aggregate_how = "magnitude")

write.csv(liana_res, "results/liana_consensus_ccc.csv", row.names = FALSE)
