## ArchR 单细胞 ATAC-seq Arrow 文件创建、TSS 富集与片段数质控
library(ArchR)

addArchRGenome("hg38")
addArchRThreads(threads = 4)

# 每个样本一个 fragments.tsv.gz；命名即 Arrow 文件与 Sample 标签
inputFiles <- c("Data/fragments_sample1.tsv.gz", "Data/fragments_sample2.tsv.gz")
names(inputFiles) <- c("sample1", "sample2")

# minTSS/minFrags 在建 Arrow 时即完成细胞级质控，并预存 Tile 与 GeneScore 矩阵
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 1,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchRProject(ArrowFiles, copyArrows = FALSE)
proj@cellColData$CellID <- rownames(proj@cellColData)
if (!("Batch" %in% colnames(proj@cellColData))) {
  proj@cellColData$Batch <- proj@cellColData$Sample
}

saveRDS(proj, "results/archr_project.rds")
