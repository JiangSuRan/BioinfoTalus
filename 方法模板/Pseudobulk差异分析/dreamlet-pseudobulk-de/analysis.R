## dreamlet 单细胞样本级 Pseudobulk 聚合与跨亚群差异表达建模
library(Seurat)
library(SingleCellExperiment)
library(dreamlet)
library(dplyr)

pbmc <- readRDS("Data/seurat_inh.rds")
metadata <- read.csv("Data/donor_metadata.csv")

# 将供体级临床表按 projid 关联到细胞级 metadata，供 pseudobulk 设计与协变量使用
data <- pbmc@meta.data[, "projid", drop = FALSE]
data2 <- left_join(data, metadata, by = "projid")
row.names(data2) <- row.names(data)
pbmc <- AddMetaData(pbmc, data2)

DefaultAssay(pbmc) <- "RNA"
pbmc.sce <- as.SingleCellExperiment(pbmc)

# 低表达基因在样本级 pseudobulk 中几乎不贡献信息且膨胀多重检验负担
pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]

# 按亚群 × 供体聚合为样本级计数，避免把同一供体的细胞当作独立重复
pb <- aggregateToPseudoBulk(
  pbmc.sce,
  assay = "counts",
  cluster_id = "cell_type_high_resolution",
  sample_id = "projid",
  verbose = FALSE
)

# 比较设计：test 为来源研究设定的脆弱抑制性亚群，baseline 为其余抑制性亚群
lst <- list(
  test = c("Inh CUX2 MSR1", "Inh L3-5 SST MAFB", "Inh LAMP5 RELN", "Inh L1 PAX6 CA4", "Inh PTPRK FAM19A1"),
  baseline = c("Inh PVALB CA8 (Chandelier)", "Inh RYR3 TSHZ2",
               "Inh PVALB HTR4", "Inh ENOX2 SPHKAP",
               "Inh VIP CLSTN2", "Inh GPC5 RIT2",
               "Inh VIP ABI3BP", "Inh PVALB SULF1",
               "Inh ALCAM TRPM3", "Inh L5-6 SST TH",
               "Inh LAMP5 NRG1 (Rosehip)", "Inh L1-6 LAMP5 CA13", "Inh L6 SST NPY",
               "Inh VIP TSHZ2", "Inh SORCS1 TTN", "Inh FBN2 EPB41L4A",
               "Inh L1-2 PAX6 SCGN", "Inh L5-6 PVALB STON2", "Inh VIP THSD7B", "Inh SGCD PDE3A")
)

fit <- dreamletCompareClusters(pb, lst, method = "fixed")

res <- topTable(fit, coef = "compare", number = Inf)

write.csv(res, "results/dreamlet_pseudobulk_de.csv")
