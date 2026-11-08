## ktplots CellPhoneDB 配受体互作显著性点图
library(ktplots)
library(ggplot2)

means <- read.delim("Data/cellphonedb_means.txt", check.names = FALSE)
pvals <- read.delim("Data/cellphonedb_pvalues.txt", check.names = FALSE)

# 目标细胞对列子集（发送方|接收方命名对齐）
means_Cancer_Treg <- means[, grep(".*MP.*Treg", colnames(means))]
pvals_Cancer_Treg <- pvals[, colnames(means_Cancer_Treg)]
means_Cancer_Treg <- cbind(pvals[, 1:11], means_Cancer_Treg)
pvals_Cancer_Treg <- cbind(pvals[, 1:11], pvals_Cancer_Treg)

# plot_cpdb：点大小 = 均值，红色边框 = 显著
p <- plot_cpdb(
  cell_type1 = "MP1_On",
  cell_type2 = "Treg",
  scdata = NULL,
  ident_key = "celltype",
  means = means_Cancer_Treg,
  pvals = pvals_Cancer_Treg,
  split.by = "orig.ident",
  cluster_rows = FALSE,
  scale = FALSE,
  keep_significant_only = TRUE,
  genes = NULL,
)

ggsave(filename = "results/ktplots_cpdb_dotplot.pdf", plot = p, width = 8, height = 10, units = "in")
