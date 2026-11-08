## AUCell 单细胞基因集排名构建、AUC 活性打分与阈值划分
library(AUCell)
library(SummarizedExperiment)
library(dplyr)

sce <- readRDS("Data/blood_tcells_sce.rds")
sig_df <- read.delim("Data/signature_genes.tsv", header = FALSE, stringsAsFactors = FALSE)

geneSets <- list(signature = sig_df$V1)

# 排名基于标准化表达而非原始计数，避免文库深度主导基因排序
cells_rankings <- AUCell_buildRankings(assay(sce, "norm_exprs"))
# aucMaxRank 取前 5% 基因，使打分只对高表达尾部敏感，对应基因开/关式激活假设
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist = FALSE, nCores = 1, assign = TRUE)

sce$AUC_signature <- SummarizedExperiment::assay(cells_AUC, "AUC")[1, ]
sce$AUC_signature.pred <- ifelse(colnames(sce) %in% cells_assignment$signature$assignment, "Positive", "Other")

# 患者与组织部位水平汇总阳性细胞频率，供组间细胞状态比例比较
info <- as.data.frame(colData(sce))
freq <- info %>%
  group_by(patient.uid, loc) %>%
  summarise(
    patient = first(patient),
    dataset = first(dataset),
    cancerType = first(cancerType),
    freq = sum(AUC_signature.pred == "Positive") / n(),
    N.all = n(),
    .groups = "drop"
  )

write.table(freq, file = "results/aucell_patient_freq.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(sce, "results/sce_with_aucell.rds")
