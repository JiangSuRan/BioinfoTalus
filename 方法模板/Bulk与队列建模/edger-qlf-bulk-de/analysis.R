## Bulk RNA-seq 原始计数矩阵 edgeR TMM 归一化与准似然 F 检验差异分析
library(edgeR)
library(limma)

counts_mat <- as.matrix(readRDS("Data/bulk_raw_counts.rds"))
meta <- read.csv("Data/bulk_metadata.csv", stringsAsFactors = FALSE)

# 严格对齐计数矩阵样本列与分组元数据顺序
common_samples <- intersect(colnames(counts_mat), meta$Sample)
counts_mat <- counts_mat[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$Sample), , drop = FALSE]
group <- factor(meta$Group, levels = c("Control", "Treat"))

# 过滤低表达噪音基因并执行 TMM 文库大小与组成偏倚校正
y <- DGEList(counts = counts_mat, group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 准似然框架（glmQLFit/glmQLFTest）对离散度不确定性建模更保守，优于 LRT
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
cont <- makeContrasts(Treat - Control, levels = design)
qlf <- glmQLFTest(fit, contrast = cont)

deg_tab <- topTags(qlf, n = Inf, adjust.method = "BH")$table
deg_tab$gene <- rownames(deg_tab)
write.csv(deg_tab, "results/edger_qlf_deg.csv", row.names = FALSE)
