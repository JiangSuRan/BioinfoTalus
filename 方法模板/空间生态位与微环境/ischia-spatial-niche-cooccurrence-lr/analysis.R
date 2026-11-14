## ISCHIA 空间组成生态位（Composition Cluster）聚类、细胞共现与富集配体-受体分析
library(Seurat)
library(ISCHIA)
library(dplyr)

visium_obj <- readRDS("Data/visium_deconvolved.rds")
deconv_mat <- readRDS("Data/deconv_weights_matrix.rds")
lr_pairs <- read.csv("Data/lr_pairs_human.csv")

# 基于反卷积细胞比例矩阵对空间 Spot 进行组成聚类，划定具有相似细胞配比的生态位微环境（CC）
visium_obj <- Composition.cluster(visium_obj, deconv_mat, k = 6)

enriched_cc4 <- Composition_cluster_enrichedCelltypes(visium_obj, "CC4", deconv_mat)
cooccur_res <- spatial.celltype.cooccurence(
  visium_obj,
  deconv.prob.mat = deconv_mat,
  COOCCUR.Threshold = 0.2,
  Condition = "Treat"
)

all_lr_genes <- intersect(unique(c(lr_pairs$ligand, lr_pairs$receptor)), rownames(visium_obj))
enriched_lrs <- Enriched.LRs(
  object = visium_obj,
  cc.sel = c("CC4"),
  idents = unique(visium_obj$orig.ident),
  all.LR.genes = all_lr_genes,
  LR.pairs = lr_pairs,
  pval.thresh = 0.05,
  logfc.thresh = 0.25
)

saveRDS(visium_obj, "results/ischia_niche_seurat.rds")
write.csv(enriched_cc4, "results/ischia_cc4_enriched_celltypes.csv", row.names = FALSE)
write.csv(enriched_lrs, "results/ischia_enriched_lrs.csv", row.names = FALSE)
