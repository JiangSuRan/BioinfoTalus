## 差异基因 GO 与 KEGG 富集分析
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggnewscale)
library(ggplot2)
library(pathview)

# 以差异基因为例进行富集分析
load("Data/DEG_results.Rda")
# 筛选差异显著的基因
DEG <- DEG %>% filter(!is.na(padj), padj < 0.05, log2FoldChange > 2)
# 沿用正向 log2FC > 2；方向以 DEG_results.Rda 保存的 contrast 为准

# SYMBOL 由差异分析中的 gene_id 注释提供，不把 Ensembl 行名当 SYMBOL
if (!all(c("gene_id", "SYMBOL") %in% names(DEG))) {
  stop("请先运行更新后的差异分析模板，提供 gene_id 和 SYMBOL 列")
}
genelist <- bitr(unique(DEG$SYMBOL), fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)
DEG_mapped <- inner_join(DEG, genelist, by = "SYMBOL")
gene_entrez <- unique(DEG_mapped$ENTREZID)
message("输入基因: ", n_distinct(DEG$gene_id), "；映射成功基因: ",
        n_distinct(DEG_mapped$gene_id), "；ENTREZ: ", length(gene_entrez))
if (length(gene_entrez) == 0) stop("没有可用于富集的 ENTREZID")

# 沿用原模板的数据库默认背景；具体研究需确认背景后再运行 ORA
# GO 富集阈值（BH 校正下的 p 值与 q 值）
p_cutoff <- 0.05
q_cutoff <- 0.05

# GO 三本体柱状图渐变配色（浅→深，映射 -log10(p.adjust)）
enrich_bp <- c("#DCEAF7", "#6BAED6", "#08519C")  # BP：Blues 系
enrich_cc <- c("#E5F5E0", "#74C476", "#238B45")  # CC：Greens 系
enrich_mf <- c("#FEE6CE", "#FDAE6B", "#E6550D")  # MF：Oranges 系

# GO：ont="all" 一次跑完 BP/CC/MF（内部按本体分别检验与 BH 校正），
# 单一 ego 对象同时供下游 cnetplot 与分本体柱状图使用
ego <- enrichGO(gene = gene_entrez,
                OrgDb = org.Hs.eg.db,
                ont = "all",
                pAdjustMethod = "BH",
                pvalueCutoff = p_cutoff,
                qvalueCutoff = q_cutoff,
                readable = TRUE)

ego_res <- as.data.frame(ego)
save(ego,ego_res,file = "Data/GO_results.Rda")


# KEGG
# 需要联网
options(timeout = 300)
kk <- enrichKEGG(gene = gene_entrez,
                 organism = 'hsa', # 物种
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.1)
kk_res <- as.data.frame(kk)
save(kk,kk_res,file = "Data/KEGG_results.Rda")

#网络图
# 生成绘图所需的foldChange List
# 没有foldChange数据就不执行如下代码，且绘图时不传入foldChange参数
# 一对多/多对一映射可能使颜色值不唯一，不能静默择优或取首条
if (anyDuplicated(DEG_mapped$ENTREZID) || anyDuplicated(DEG_mapped$SYMBOL)) {
  stop("映射含重复 ID；富集结果已保存，绘制表达着色图前需明确重复 ID 的处理规则")
}
List <- DEG_mapped$log2FoldChange
names(List) <- DEG_mapped$ENTREZID
List <- sort(List, decreasing = TRUE)

#GO
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
#KEGG
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)

# 气泡图和条形图
# GO 分本体柱状图：各本体最显著的前 10 条，颜色按 -log10(p.adjust) 渐变
go_bar_plot <- function(enrich_df, plot_title, bar_colors, top_n = 10) {
  enrich_df %>%
    arrange(p.adjust, pvalue) %>%
    slice_head(n = top_n) %>%
    mutate(
      neg_log10_padj = -log10(p.adjust),
      Description = str_wrap(Description, width = 40),        # 描述过长自动换行
      Description = fct_reorder(Description, neg_log10_padj)  # 最显著条目排顶部
    ) %>%
    ggplot(aes(x = Count, y = Description, fill = neg_log10_padj)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradientn(colors = bar_colors, name = "-log10(p.adjust)") +
    labs(title = plot_title, x = "Gene Number", y = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          axis.text.y = element_text(size = 9))
}

# 各本体分别按校正 p 值排序；显式保留图对象以免保存成最后一张图
p_bp <- go_bar_plot(filter(ego_res, ONTOLOGY == "BP"), "GO Biological Process", enrich_bp)
p_cc <- go_bar_plot(filter(ego_res, ONTOLOGY == "CC"), "GO Cellular Component", enrich_cc)
p_mf <- go_bar_plot(filter(ego_res, ONTOLOGY == "MF"), "GO Molecular Function", enrich_mf)
print(p_bp)
print(p_cc)
print(p_mf)
# 需要保存时取消注释
# ggsave("GO_BP_barplot.pdf", p_bp, width = 7, height = 5)

# KEGG
kk_df <- as.data.frame(kk)
dotplot(kk, showCategory=15, orderBy="GeneRatio", label_format=70)
barplot(kk, showCategory=15, drop=TRUE, label_format=70)

# KEGG通路可视化
pathway <- pathview(gene.data=unique(DEG_mapped$SYMBOL), # 基因列表
                  pathway.id="hsa04820", # 需要可视化的通路名
                  species="hsa", # 物种
                  gene.idtype="SYMBOL", # 基因ID类型 
                  kegg.native=TRUE)

# KEGG通路+差异表达可视化
deg_kegg <- DEG_mapped[, c("SYMBOL", "log2FoldChange")]
deg_kegg <- column_to_rownames(deg_kegg, "SYMBOL")
pathway<-pathview(gene.data=deg_kegg, # 基因列表
                  pathway.id="hsa04820", # 通路名
                  species="hsa", # 物种
                  gene.idtype="SYMBOL", # 基因ID类型 
                  kegg.native=TRUE)

