## chromVAR motif 偏差 Z-score（GC 偏倚校正 + motifmatchr 匹配 + 偏差变异度）
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(chromVARmotifs)

peaks <- getPeaks("Data/peaks.bed", sort_peaks = TRUE)
sampling_depth <- as.numeric(read.delim("Data/sampling_depth.txt", header = FALSE))
my_counts_matrix <- as.matrix(read.delim("Data/counts.txt", header = FALSE))
rownames(my_counts_matrix) <- as.matrix(read.delim("Data/peak_names.txt", header = FALSE))
colnames(my_counts_matrix) <- as.matrix(read.delim("Data/cell_names.txt", header = FALSE))

fragment_counts <- SummarizedExperiment(
  assays = list(counts = my_counts_matrix),
  rowRanges = peaks,
  colData = sampling_depth
)
colnames(colData(fragment_counts)) <- "depth"

# GC 含量偏倚校正（motif 偏差计算的前置）
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)

# motif 匹配（chromVARmotifs human_pwms_v2）
data("human_pwms_v2")
motif_ix <- matchMotifs(human_pwms_v2, fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)

# 偏差计算：每细胞每 motif 的校正后可及性偏差 Z-score 与变异度
dev <- computeDeviations(object = fragment_counts, annotations = motif_ix)
variability <- computeVariability(dev)

write.csv(deviations(dev), "results/chromvar_deviations.csv")
write.csv(deviationScores(dev), "results/chromvar_deviationScores.csv")
write.csv(variability, "results/chromvar_variability.csv")
