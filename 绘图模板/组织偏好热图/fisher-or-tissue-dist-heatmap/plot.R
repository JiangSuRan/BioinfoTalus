## Fisher 精确检验组织富集 Odds Ratio 热图
library(RColorBrewer)
library(sscVis)

meta <- readRDS("Data/cell_tissue_metadata.rds")

tab <- table(as.character(meta$cell_type), as.character(meta$tissue))
tot <- sum(tab)
rsum <- rowSums(tab)
csum <- colSums(tab)
or_mat <- matrix(NA_real_, nrow = nrow(tab), ncol = ncol(tab), dimnames = dimnames(tab))

for (i in seq_len(nrow(tab))) {
  for (j in seq_len(ncol(tab))) {
    m <- matrix(c(tab[i, j], rsum[i] - tab[i, j], csum[j] - tab[i, j], tot - rsum[i] - csum[j] + tab[i, j]), nrow = 2)
    or_mat[i, j] <- fisher.test(m)$estimate
  }
}

# 将极端优势比截断在 2，避免少数高富集亚群压缩整体色阶动态范围
or_capped <- pmin(round(or_mat, 2), 2)

sscVis::plotMatrix.simple(
  or_capped,
  out.prefix = "results/tissue_or_dist",
  show.number = TRUE,
  do.clust = FALSE,
  clust.column = FALSE,
  clust.row = FALSE,
  show.dendrogram = TRUE,
  waterfall.row = FALSE,
  exp.name = expression(italic(OR)),
  z.hi = 2,
  palatte = brewer.pal(9, "YlGnBu"),
  pdf.width = 6,
  pdf.height = 8
)
