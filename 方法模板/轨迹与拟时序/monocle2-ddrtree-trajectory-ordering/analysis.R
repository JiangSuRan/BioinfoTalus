## Monocle 2 DDRTree 拟时序降维、状态排序与分支差异 BEAM 分析
library(Seurat)
library(monocle)
library(dplyr)

sc_obj <- readRDS("Data/seurat_subset.rds")

# Seurat v5 显式通过 LayerData 读取稀疏 counts 并同步构建 CellDataSet 所需的注释框
expr_mat <- LayerData(sc_obj, assay = "RNA", layer = "counts")
pd <- new("AnnotatedDataFrame", data = sc_obj@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat)))

cds <- newCellDataSet(
  expr_mat,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff_res <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~celltype")
ordering_genes <- row.names(subset(diff_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)

# 针对分支点 1 运行 BEAM 检验以识别决定细胞命运分叉的分支依赖基因
beam_res <- BEAM(cds, branch_point = 1, cores = 1)
beam_res <- beam_res %>% arrange(qval)

saveRDS(cds, "results/monocle2_cds.rds")
write.csv(beam_res, "results/monocle2_beam_genes.csv", row.names = TRUE)
