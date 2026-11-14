## Beyondcell 单细胞药物敏感性特征打分、技术混杂回归与药物排名
library(beyondcell)
library(Seurat)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- ScaleData(sc_obj, verbose = FALSE)
# beyondcell 当前读取经典 Assay 槽位，若传入 Seurat v5 Assay5 先转回经典 Assay 结构
sc_obj[["RNA"]] <- as(object = sc_obj[["RNA"]], Class = "Assay")

drug_gs <- GetCollection(PSc)
bc <- bcScore(sc_obj, drug_gs, expr.thres = 0.1)
bc <- bcUMAP(bc, pc = 10, k.neighbors = 4, res = 0.2)
# 剥除单细胞检测基因数对药敏得分空间的驱动效应后重算药敏空间 UMAP
bc <- bcRegressOut(bc, vars.to.regress = "nFeature_RNA")
bc <- bcUMAP(bc, pc = 10, k.neighbors = 20, res = 0.2)

bc <- bcRanks(bc, idents = "condition")
bc <- bcRanks(bc, idents = "bc_clusters_res.0.2", extended = FALSE)

saveRDS(bc, "results/beyondcell_scored.rds")
write.table(bc@ranks$condition, "results/beyondcell_condition_ranks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
