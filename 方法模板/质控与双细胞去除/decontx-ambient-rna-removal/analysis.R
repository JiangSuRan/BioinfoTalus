## DecontX 环境游离 RNA 污染估计、矫正表达写回与低污染细胞过滤
library(decontX)
library(Seurat)

sc_obj <- readRDS("Data/seurat_obj.rds")
counts_mat <- LayerData(sc_obj, assay = "RNA", layer = "counts")

decontx_res <- decontX(counts_mat)
sc_obj$Contamination <- decontx_res$contamination
# 同时保存 DecontX 估计的去污染计数层供下游参考，避免直接覆写原始 counts
sc_obj[["decontXcounts"]] <- CreateAssayObject(counts = round(decontx_res$decontXcounts))

sc_filtered <- subset(sc_obj, subset = Contamination < 0.2)
saveRDS(sc_filtered, "results/seurat_decontx_filtered.rds")
