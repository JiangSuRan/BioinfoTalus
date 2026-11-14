## WGCNA 共表达聚类树模块色带图与模块-临床性状相关性标注热图
library(WGCNA)

wgcna_res <- readRDS("results/wgcna_modules.rds")

# 基因系统发育树与动态剪切/合并前后模块色带对照图
pdf("results/wgcna_dendro_modules.pdf", width = 10, height = 5.5)
plotDendroAndColors(
  dendro = wgcna_res$geneTree,
  colors = cbind(wgcna_res$dynamicColors, wgcna_res$moduleColors),
  groupLabels = c("Dynamic Cut", "Merged Modules"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram and Co-expression Module Assignment"
)
dev.off()

cor_mat <- wgcna_res$moduleTraitCor
pval_mat <- wgcna_res$moduleTraitPvalue
text_matrix <- paste0(signif(cor_mat, 2), "\n(", signif(pval_mat, 1), ")")
dim(text_matrix) <- dim(cor_mat)

# yLabels 保留 "ME" 前缀启用色块注释，细胞内同时标注相关系数与显著性 p 值
pdf("results/wgcna_module_trait_heatmap.pdf", width = 8.5, height = 6.5)
par(mar = c(6, 9, 3, 2))
labeledHeatmap(
  Matrix = cor_mat,
  xLabels = colnames(cor_mat),
  yLabels = rownames(cor_mat),
  ySymbols = rownames(cor_mat),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.72,
  cex.lab = 0.85,
  zlim = c(-1, 1),
  main = "Module-Trait Pearson Correlations"
)
dev.off()
