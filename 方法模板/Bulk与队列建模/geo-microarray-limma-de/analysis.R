## GEO / 芯片表达矩阵探针映射、分位数归一化与 limma 经验贝叶斯差异分析
library(limma)
library(dplyr)

expr_probe <- as.matrix(readRDS("Data/microarray_probe_expr.rds"))
probe_map <- read.table("Data/probe_annotation.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
meta <- read.table("Data/sample_metadata.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 剔除缺失或多基因共注释（///）探针，再用 avereps 按基因均值合并多探针
probe_clean <- probe_map %>%
  filter(!is.na(SYMBOL), SYMBOL != "", !grepl("///", SYMBOL)) %>%
  mutate(ID = as.character(ID), SYMBOL = trimws(SYMBOL))
common_probes <- intersect(rownames(expr_probe), probe_clean$ID)
expr_sub <- expr_probe[common_probes, , drop = FALSE]
probe_symbols <- probe_clean$SYMBOL[match(common_probes, probe_clean$ID)]
expr_gene <- avereps(expr_sub, ID = probe_symbols)

# 对齐样本分组与表达矩阵列顺序，分位数归一化消除芯片阵列间技术偏移
common_samples <- intersect(colnames(expr_gene), meta$Sample)
meta <- meta[match(common_samples, meta$Sample), , drop = FALSE]
expr_norm <- normalizeBetweenArrays(expr_gene[, common_samples, drop = FALSE], method = "quantile")

group <- factor(meta$Group, levels = c("Control", "Tumor"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# trend = TRUE 与 robust = TRUE 兼顾平均信号强度方差趋势及离群样本稳健性
fit <- lmFit(expr_norm, design)
cont <- makeContrasts(Tumor - Control, levels = design)
fit2 <- contrasts.fit(fit, cont) %>%
  eBayes(trend = TRUE, robust = TRUE)

deg_tab <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
deg_tab$gene <- rownames(deg_tab)
write.csv(deg_tab, "results/limma_microarray_deg.csv", row.names = FALSE)
