## GeneNMF 跨样本一致表达元程序发现、Hallmark 功能注释与 UCell 单细胞打分
library(Seurat)
library(GeneNMF)
library(fgsea)
library(UCell)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
DefaultAssay(sc_obj) <- "RNA"
seu_list <- SplitObject(sc_obj, split.by = "sample")

# 在各样本内跨多个分解秩运行 NMF，减少单样本自由度与批次扰动的影响
nmf_programs <- multiNMF(seu_list, assay = "RNA", k = 4:9, nfeatures = 1000)

meta_programs <- getMetaPrograms(
  nmf_programs,
  metric = "cosine",
  weight.explained = 0.5,
  nMP = 6,
  min.confidence = 0.7
)

mp_gsea <- lapply(meta_programs$metaprograms.genes, function(genes) {
  runGSEA(genes, universe = rownames(sc_obj), category = "H")
})

mp_genes <- meta_programs$metaprograms.genes
sc_obj <- AddModuleScore_UCell(sc_obj, features = mp_genes, ncores = 1, name = "")

saveRDS(meta_programs, "results/genenmf_metaprograms.rds")
saveRDS(mp_gsea, "results/genenmf_gsea_hallmark.rds")
saveRDS(sc_obj, "results/seurat_with_mp_ucell.rds")
