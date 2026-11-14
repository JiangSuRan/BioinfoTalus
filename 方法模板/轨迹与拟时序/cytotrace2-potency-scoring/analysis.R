## CytoTRACE2 单细胞分化潜能预测与 Seurat 元数据回写
library(Seurat)
library(CytoTRACE2)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
counts_mat <- GetAssayData(sc_obj, assay = "RNA", layer = "counts")

# 使用原始 counts 与显式随机种子预测分化潜能，供下游轨迹根节点先验确定
ct2_res <- cytotrace2(
  counts_mat,
  is_seurat = FALSE,
  slot_type = "counts",
  species = "human",
  parallelize_models = FALSE,
  parallelize_smoothing = FALSE,
  seed = 1234
)

sc_obj <- AddMetaData(sc_obj, metadata = ct2_res)
saveRDS(ct2_res, "results/cytotrace2_result.rds")
saveRDS(sc_obj, "results/seurat_with_cytotrace2.rds")
