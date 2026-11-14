## TCGA 与 GTEx FPKM 对数尺度统一、GSVA 2.x 通路评分及 limma 差异检验
library(limma)
library(GSVA)
library(BiocParallel)
library(clusterProfiler)

tcga_expr <- as.matrix(readRDS("Data/tcga_log2fpkm1.rds"))
gtex_expr <- as.matrix(readRDS("Data/gtex_log2fpkm001.rds"))

# Xena Toil 平台中 GTEx 为 log2(FPKM+0.001) 而 TCGA 为 log2(FPKM+1)，必须还原统一偏移后再合并
gtex_fpkm <- pmax(2^gtex_expr - 0.001, 0)
gtex_aligned <- log2(gtex_fpkm + 1)

common_genes <- intersect(rownames(tcga_expr), rownames(gtex_aligned))
merged_expr <- cbind(
  tcga_expr[common_genes, , drop = FALSE],
  gtex_aligned[common_genes, , drop = FALSE]
)
expr_norm <- normalizeBetweenArrays(merged_expr, method = "quantile")

group <- factor(
  c(rep("Tumor", ncol(tcga_expr)), rep("Normal", ncol(gtex_aligned))),
  levels = c("Normal", "Tumor")
)

gmt_df <- read.gmt("Data/hallmark_genesets.gmt")
pathway_list <- split(gmt_df$gene, gmt_df$term)

# 使用 GSVA 2.x 参数对象接口评估各样本连续表达通路活性
gsva_param <- gsvaParam(
  exprData = expr_norm,
  geneSets = pathway_list,
  minSize = 5,
  maxSize = 500,
  kcdf = "Gaussian",
  maxDiff = TRUE
)
gsva_mat <- gsva(gsva_param, verbose = FALSE, BPPARAM = SerialParam(progressbar = FALSE))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(gsva_mat, design)
cont <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont) %>%
  eBayes()

pathway_diff <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
pathway_diff$pathway <- rownames(pathway_diff)

saveRDS(gsva_mat, "results/tcga_gtex_gsva_scores.rds")
write.csv(pathway_diff, "results/tcga_gtex_gsva_limma.csv", row.names = FALSE)
