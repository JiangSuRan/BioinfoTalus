## 单细胞 Seurat 质控指标计算、DoubletFinder 双细胞识别与 MAD 异常值过滤
library(Seurat)
library(DoubletFinder)
library(scater)
library(stringr)

# 线粒体、核糖体、血红蛋白与血小板基因占比是区分低质量细胞和双细胞的通用质控协变量
counts <- Read10X("Data/filtered_feature_bc_matrix")
seu <- CreateSeuratObject(counts)
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent.mt")
seu <- PercentageFeatureSet(seu, "^RP[SL]", col.name = "percent.ribo")
seu <- PercentageFeatureSet(seu, "^HB[^(P)]", col.name = "percent.hb")
seu <- PercentageFeatureSet(seu, "PECAM1|PF4", col.name = "percent.plat")

# 初步过滤去除基因数过少与线粒体/核糖体占比过高的细胞，避免极端细胞干扰降维与 pK 扫描
seu <- subset(seu, subset = nFeature_RNA > 500 & percent.mt < 1 & percent.ribo < 1)

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, dims = 1:15)

# DoubletFinder 需先在相同 PCs 上做 pN-pK 参数扫描，BCmetric 峰值处即为该数据最佳 pK
sweep.res.list <- paramSweep(seu, PCs = 1:15, num.cores = 5)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# 预期双细胞数按 10x 试剂 multiplet 率随细胞数外推，再扣除同源双细胞比例
annotations <- seu$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(1/100 * (ncol(seu)/1000) * ncol(seu))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

seu <- doubletFinder(
  seu,
  PCs = 1:15,
  pN = 0.25,
  pK = mpK,
  nExp = nExp_poi.adj,
  reuse.pANN = FALSE
)

# 对数尺度上 3 MAD 稳健过滤，阈值随样本自身分布自适应，避免跨样本硬阈值漂移
seu.nCount <- isOutlier(seu$nCount_RNA, nmads = 3, type = "lower", log = TRUE)
seu.nFeatures <- isOutlier(seu$nFeature_RNA, nmads = 3, type = "lower", log = TRUE)
seu.mito <- isOutlier(seu$percent.mt, nmads = 3, type = "higher", log = TRUE)
seu.doublet <- seu@meta.data[, str_detect(colnames(seu@meta.data), "DF.classifications")] == "Doublet"
cells.out <- seu.mito | seu.nFeatures | seu.nCount | seu.doublet

seu <- seu[, !cells.out]

saveRDS(seu, "results/seurat_qc_filtered.rds")
