## TCGAbiolinks GDC 队列获取、蛋白编码 TPM 矩阵整理与去重
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tibble)

# GDC 查询与下载（STAR - Counts 工作流）
query <- GDCquery(
  project = "TCGA-UCEC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api", files.per.chunk = 100)
exprs.sce <- GDCprepare(query)
save(exprs.sce, file = "results/TCGA-UCEC_mRNA.Rdata")

# 仅蛋白编码基因 × 原发肿瘤样本；提取 unstranded 计数与 TPM
se_mrna <- exprs.sce[rowData(exprs.sce)$gene_type == "protein_coding",
                     colData(exprs.sce)$sample_type == "Primary Tumor"]
exprs_mrna_count <- assay(se_mrna, "unstranded")
exprs_mrna_tpm <- assay(se_mrna, "tpm_unstrand")

sample_id <- colData(se_mrna)$patient
exprs_mrna_count <- cbind(data.frame(sample_id), as.data.frame(t(exprs_mrna_count)))
exprs_mrna_tpm <- cbind(data.frame(sample_id), as.data.frame(t(exprs_mrna_tpm)))

# 患者多播片去重（保留均值最高条码）
exprs_mrna_count <- exprs_mrna_count %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "sample_id") %>%
  as.data.frame()

exprs_mrna_tpm <- exprs_mrna_tpm %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "sample_id") %>%
  as.data.frame()

write.csv(exprs_mrna_count, "results/tcga_ucec_counts.csv")
write.csv(exprs_mrna_tpm, "results/tcga_ucec_tpm.csv")
