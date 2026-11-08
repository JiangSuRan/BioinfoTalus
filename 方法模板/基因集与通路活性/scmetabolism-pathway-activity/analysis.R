## scMetabolism 单细胞 KEGG 代谢通路活性量化（AUCell 模式）
library(Seurat)
library(scMetabolism)

objN <- readRDS("Data/panB_scRNA_processed_data.rds")

# 限定目标 B 细胞亚群，使通路活性在同源细胞类型内解释
Idents(objN) <- "celltype"
objectTP <- subset(objN, idents = c("B.08.ITGB1+SwBm", "B.09.DUSP4+AtM"))

# AUCell 模式逐细胞做全基因组排序，计算开销大，抽样 10000 细胞以控制运行时间与内存；固定种子保证抽样可重复
set.seed(100)
cellid <- sample(1:ncol(objectTP), 10000, replace = FALSE)
objectTPN <- objectTP[, cellid]

# imputation = FALSE 避免插值引入虚假平滑，ncores 按本机核数设定
countexp.Seurat <- sc.metabolism.Seurat(
  obj = objectTPN,
  method = "AUCell",
  imputation = FALSE,
  ncores = 4,
  metabolism.type = "KEGG"
)

saveRDS(countexp.Seurat, file = "results/seurat_scmetabolism_kegg.rds")
