## Pseudobulk 样本×亚群平均表达三维 PCA 散点图（scatterplot3d）
library(dplyr)
library(forcats)
library(scatterplot3d)

# 依据 scatterplot3d 的坐标比例尺与倾角投影 xy/xz/yz 背景网格
addgrids3d <- function(x, y = NULL, z = NULL, grid = TRUE,
                       col.grid = "grey", lty.grid = par("lty"),
                       lab = par("lab"), lab.z = mean(lab[1:2]),
                       scale.y = 1, angle = 40,
                       xlim = NULL, ylim = NULL, zlim = NULL) {
  if (inherits(x, c("matrix", "data.frame"))) {
    x <- as.data.frame(x)
    y <- unlist(x[, 2])
    z <- unlist(x[, 3])
    x <- unlist(x[, 1])
  }
  p.lab <- par("lab")
  angle_norm <- (angle %% 360) / 90
  yz.f <- scale.y * abs(if (angle_norm < 1) angle_norm else if (angle_norm > 3) angle_norm - 4 else 2 - angle_norm)
  yx.f <- scale.y * (if (angle_norm < 2) 1 - angle_norm else angle_norm - 3)

  x.range <- range(x[is.finite(x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 * lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x / x.scal
  x.max <- ceiling(range(x.prty)[2] / x.scal)
  x.min <- floor(range(x.prty)[1] / x.scal)
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2] / x.scal))
    x.min <- min(x.min, floor(xlim[1] / x.scal))
  }

  y.range <- range(y[is.finite(y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 * lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add) / y.scal
  y.max <- (max(y.prty) - y.add) / y.scal
  if (!is.null(ylim)) {
    y.max <- max(y.max, ceiling((ylim[2] - y.add) / y.scal))
  }

  z.range <- range(z[is.finite(z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 * lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z / z.scal
  z.max <- ceiling(range(z.prty)[2] / z.scal)
  z.min <- floor(range(z.prty)[1] / z.scal)
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2] / z.scal))
    z.min <- min(z.min, floor(zlim[1] / z.scal))
  }

  if ("xy" %in% grid || identical(grid, TRUE)) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, i + (yx.f * y.max), yz.f * y.max + z.max, col = col.grid, lty = lty.grid)
    i <- z.min:z.max
    segments(x.min + (yx.f * y.max), (yz.f * y.max) + i, x.max + (yx.f * y.max), (yz.f * y.max) + i, col = col.grid, lty = lty.grid)
  }
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.min + (i * yx.f), i * yz.f + z.max, col = col.grid, lty = lty.grid)
    i <- z.min:z.max
    segments(x.min + (yx.f * y.max), (yz.f * y.max) + i, x.min, i, col = col.grid, lty = lty.grid)
  }
}

avg_expr <- readRDS("Data/pseudobulk_avg_expr.rds")
pca_meta <- readRDS("Data/pseudobulk_pca_meta.rds")

pca <- prcomp(t(avg_expr))
pca_df <- as.data.frame(pca$x[, 1:3, drop = FALSE])
pca_df$id <- rownames(pca_df)
pca_df <- pca_df %>%
  left_join(pca_meta, by = "id") %>%
  mutate(
    group = as_factor(group),
    cell_type = as_factor(cell_type)
  )

group_levels <- levels(pca_df$group)
celltype_levels <- levels(pca_df$cell_type)
shape_palette <- setNames(rep_len(c(11, 16, 17, 15, 18, 8), length(group_levels)), group_levels)
color_palette <- setNames(
  rep_len(
    c("#1F78B4", "#E31A1C", "#FF7F00", "#33A02C", "#6A3D9A", "#A6CEE3", "#FB9A99", "#FDBF6F", "#B2DF8A", "#CAB2D6", "#B15928"),
    length(celltype_levels)
  ),
  celltype_levels
)

pdf("results/pseudobulk_3d_pca.pdf", width = 6.5, height = 6)
scatterplot3d(
  pca_df[, c("PC1", "PC2", "PC3")],
  pch = shape_palette[as.character(pca_df$group)],
  color = color_palette[as.character(pca_df$cell_type)],
  cex.symbols = 1,
  grid = FALSE,
  box = FALSE,
  main = "3D PCA plot"
)
legend("topright", legend = celltype_levels, fill = color_palette[celltype_levels], box.col = NA)
legend("topleft", legend = group_levels, pch = shape_palette[group_levels], box.col = NA)
addgrids3d(pca_df[, c("PC1", "PC2", "PC3")], grid = c("xy", "xz", "yz"))
dev.off()
