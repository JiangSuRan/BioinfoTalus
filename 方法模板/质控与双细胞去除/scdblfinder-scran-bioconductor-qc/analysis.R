## scDblFinder 逐样本双细胞检测与 SingCellExperiment 体系质控
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(scran)

# 多样本合并对象（min.cells=3 已过滤低表达基因）
seu <- readRDS("Data/kpc_timecourse_seurat.rds")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu <- subset(seu, subset = percent.mt < 25 & nCount_RNA > 1000 & nFeature_RNA > 200)

# 逐样本 scDblFinder：dbr 期望双细胞率 0.07（10x 常规载量经验值）
sample_ids <- unique(seu$orig.ident)
dbl_colData <- NULL
for (i in sample_ids) {
  sub <- subset(seu, subset = orig.ident == i)
  sceDblF <- scDblFinder(sub@assays$RNA@counts, dbr = 0.07)
  dbl_colData <- rbind(dbl_colData, as.data.frame(colData(sceDblF)))
}

seu$is.doublet <- dbl_colData[colnames(seu), "scDblFinder.class"]
seu$scDblFinder.score <- dbl_colData[colnames(seu), "scDblFinder.score"]
seu <- subset(seu, subset = is.doublet == "singlet")

# scran 去卷积大小因子（counts 层）
sce <- as.SingleCellExperiment(seu)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
summary(sizeFactors(sce))

saveRDS(seu, "results/kpc_qc_dblfiltered.rds")
write.csv(as.data.frame(dbl_colData)[, c("scDblFinder.class", "scDblFinder.score", "scDblFinder.doublet")],
          "results/scdblfinder_results.csv")
