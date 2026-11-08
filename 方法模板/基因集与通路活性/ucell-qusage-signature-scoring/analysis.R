## UCell Mann-Whitney U 稳健基因集打分与 qusage GMT 读入
library(UCell)
library(qusage)
library(Seurat)
library(dplyr)

seu <- readRDS("Data/neutrophil_seurat.rds")
DefaultAssay(seu) <- "RNA"

# GMT 读入（qusage 读取器，兼容 MSigDB/自建签名）
signature1 <- qusage::read.gmt("Data/NI_signature_NEU_220504.gmt")

# 行名统一大写以匹配 GMT 基因大小写
input.mat <- GetAssayData(seu, assay = "RNA", layer = "data")
row.names(input.mat) <- toupper(row.names(input.mat))

# UCell：基于排名的 Mann-Whitney U 统计量，maxRank 截断降低 dropout 影响
signature.matrix.neu <- data.frame(ScoreSignatures_UCell(input.mat, features = signature1, maxRank = 3000, name = ""))

# 打分写回 metadata
seu@meta.data <- cbind(seu@meta.data, signature.matrix.neu)

write.csv(signature.matrix.neu, "results/ucell_scores.csv")
saveRDS(seu, "results/ucell_annotated.rds")
