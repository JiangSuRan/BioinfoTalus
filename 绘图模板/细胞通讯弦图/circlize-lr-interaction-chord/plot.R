## circlize 配体-受体细胞通讯方向性多层弦图（chordDiagram）
library(circlize)
library(dplyr)

lr_df <- read.csv("data/lr_links.csv", stringsAsFactors = FALSE)
cls_cols <- readRDS("data/celltype_colors.rds")

# -log10(p) 截断到 2.5 再归一化为连线宽度，未过 0.05 阈值的弱连接置 0
plot_data <- lr_df %>%
  mutate(mlogp = -log10(pval)) %>%
  mutate(mlogp = pmin(mlogp, 2.5) / 2.5)
plot_data$mlogp[plot_data$mlogp <= -log10(0.05) / 2.5] <- 0
subdata <- plot_data %>% select(from, to, mlogp)

# 扇区半径按两端扇区连线总量自适应，保证两端高度对齐
max_sec <- 10
max_sum1 <- tapply(subdata$mlogp, subdata$from, sum) %>% max(., na.rm = TRUE)
max_sum2 <- tapply(subdata$mlogp, subdata$to, sum) %>% max(., na.rm = TRUE)
if (max(c(max_sum1, max_sum2)) > 10) max_sec <- max(c(max_sum1, max_sum2))

xmax_values <- rep(max_sec, length(cls_cols)) %>% setNames(., names(cls_cols))

pdf("result/lr_chord.pdf", width = 7, height = 7)
circos.par(start.degree = 90, cell.padding = c(0, 0, 0, 0), gap.degree = 5.5)
chordDiagram(
  subdata,
  grid.col = cls_cols,
  order = names(cls_cols),
  xmax = xmax_values,
  reduce = -1,
  col = "black",
  scale = FALSE,
  annotationTrack = c("grid"),
  annotationTrackHeight = convert_height(7, "mm"),
  transparency = 0.6,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  link.arr.length = 0.05,
  link.arr.width = 0.01
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), mean(ylim), sector.name, facing = "bending.inside",
              niceFacing = TRUE, col = "white", cex = 1.6)
}, bg.border = NA)
circos.clear()
dev.off()
