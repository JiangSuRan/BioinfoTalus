## TCGA 原始计数差异分析
library(tidyverse)
library(DESeq2)
library(ggpubr)
library(ggthemes)

input_dir <- "Data/TCGA_expression_output"
counts_matrix <- readRDS(file.path(input_dir, "counts_protein_coding.rds"))
sample_meta <- read_tsv(file.path(input_dir, "sample_metadata.tsv"),
                        col_types = cols(sample_type_code = col_character()))
gene_annotation <- read_tsv(file.path(input_dir, "gene_annotation.tsv"), show_col_types = FALSE)
fdr_cutoff <- 0.05

sample_meta <- sample_meta %>% filter(keep)
if (anyDuplicated(sample_meta$sample_barcode)) stop("样本信息存在重复条码")
idx <- match(colnames(counts_matrix), sample_meta$sample_barcode)
if (anyNA(idx)) stop("矩阵样本无法在 sample_metadata 中找到")
sample_meta <- sample_meta[idx, ]
if (!all(sample_meta$sample_type_code %in% c("01", "11"))) {
  stop("此模板只对比 01 与 11；其他样本类型需明确设计后再分析")
}
group <- ifelse(sample_meta$sample_type_code == "01", "tumor", "normal")
conditions <- data.frame(group = factor(group, levels = c("tumor", "normal")),
                          row.names = sample_meta$sample_barcode)
if (any(table(conditions$group) == 0)) stop("对比需要 tumor 和 normal 两组")
if (any(!is.finite(counts_matrix)) || any(counts_matrix < 0) ||
    any(counts_matrix != floor(counts_matrix))) stop("必须使用非负整数原始计数，不用 round 修补")
if (anyDuplicated(rownames(counts_matrix)) || anyDuplicated(gene_annotation$gene_id)) {
  stop("gene_id 重复，需先核实输入")
}

#+ Check Point：原始 counts 与样本表已对齐；沿用原模板 ~ group 的独立样本设计
# 若队列含配对或协变量，应先明确设计，不能直接沿用本示例
# 保持旧版 normal/tumor 方向，正 log2FC 表示 normal 较高
contrast <- c("group", "normal", "tumor")
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = conditions,
                              design = ~ group)
dds <- DESeq(dds)
DEG_full <- as.data.frame(results(dds, contrast = contrast, alpha = fdr_cutoff))
DEG_full$gene_id <- rownames(DEG_full)
annotation_idx <- match(DEG_full$gene_id, gene_annotation$gene_id)
if (anyNA(annotation_idx)) stop("差异结果 gene_id 无法对应基因注释")
DEG_full$SYMBOL <- gene_annotation$gene_name[annotation_idx]

# 完整结果保留 NA 和未显著基因，供背景核实与火山图使用
DEG <- DEG_full %>% filter(!is.na(padj), padj < fdr_cutoff)
save(DEG, DEG_full, contrast, fdr_cutoff, file = "Data/DEG_results.Rda")
DEG_plot <- DEG_full %>% filter(!is.na(padj), !is.na(log2FoldChange))

# 火山图
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG_plot$padj < 0.05)&(DEG_plot$log2FoldChange < -logFC_cutoff)
type2 = (DEG_plot$padj < 0.05)&(DEG_plot$log2FoldChange > logFC_cutoff)
DEG_plot$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
DEG_plot$logP <- -log10(DEG_plot$padj)

# 添加基因标签信息
DEG_plot$Label = ""  # 新加一列label
# 对差异基因的p值进行从小到大的排序
DEG_plot <- DEG_plot[order(DEG_plot$padj), ]
# 将行名赋给DEG_plot$Gene
DEG_plot$Gene <- rownames(DEG_plot)
# 高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG_plot$Gene[which(DEG_plot$change == "UP")], 5)
# 低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG_plot$Gene[which(DEG_plot$change == "DOWN")], 5)
# 将up.genes和down.genes合并，并加入到Label中
DEG_plot.top5.genes <- c(as.character(up.genes), as.character(down.genes))
# 在DEG_plot$Label中为DEG_plot.top5.genes对应的基因添加标签
DEG_plot$Label[match(DEG_plot.top5.genes, DEG_plot$Gene)] <- DEG_plot.top5.genes

ggscatter(DEG_plot, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG_plot$Label,
          font.label = 8,
          repel = TRUE,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed")



