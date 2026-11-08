## DirichletReg 细胞组成回归与 MAST 零膨胀双部分差异表达
library(DirichletReg)
library(MAST)
library(Matrix.utils)
library(data.table)

seur <- readRDS("Data/uc_seurat.rds")

# 组成回归：样本 × 细胞类型计数表 + 条件协变量
counts_ct <- do.call(cbind, tapply(seur@meta.data$orig.ident, seur@meta.data$celltype, table))
covariates <- data.frame(condition = gsub("[12].*", "", rownames(counts_ct)))

counts_ct <- as.data.frame(counts_ct)
counts_ct$counts <- DR_data(counts_ct)
data <- cbind(counts_ct, covariates)
fit <- DirichReg(counts ~ condition, data)

# 提取各细胞类型 condition 效应 p 值矩阵
u <- summary(fit)
pvals <- u$coef.mat[grep("Intercept", rownames(u$coef.mat), invert = TRUE), 4]
v <- names(pvals)
pvals_mat <- matrix(pvals, ncol = length(u$varnames))
rownames(pvals_mat) <- gsub("condition", "", v[seq_len(nrow(pvals_mat))])
colnames(pvals_mat) <- u$varnames
fit$pvals <- pvals_mat

saveRDS(fit, "results/dirichlet_reg_fit.rds")
write.csv(fit$pvals, "results/dirichlet_reg_pvals.csv")

# MAST：阈值化 + zlm 零膨胀双部分模型 + 条件对比
Idents(seur) <- "celltype"
mono <- subset(seur, idents = "Monocytes")

sca <- FromMatrix(
  exprsArray = as.matrix(GetAssayData(mono, assay = "RNA", slot = "data")),
  cData = mono@meta.data
)

sca <- sca[apply(assay(sca), 1, function(x) sum(x > 0) > 10), ]
thres <- thresholdSCRNASEQExpressionClean(assay(sca), classes = rep("none", ncol(sca)), medianMarker = FALSE)
sca <- threshold(sca, cutoff = thres)

zlm_cond <- zlm(~ condition + nCount_RNA, sca)
summary_cond <- summary(zlm_cond, doLRT = "conditionUC")
mast_res <- summary_cond$datatable

fwrite(mast_res, "results/mast_zlm_condition_lrt.csv")
