## 多条件 CellChat 通讯数量/强度对比、差异网络与 rankNet 通路富集对比图
library(CellChat)
library(patchwork)
library(ComplexHeatmap)

cc_list <- readRDS("Data/cellchat_list.rds")
cellchat <- mergeCellChat(cc_list, add.names = names(cc_list), cell.prefix = TRUE)
grp_idx <- seq_along(cc_list)

p_cnt <- compareInteractions(cellchat, show.legend = FALSE, group = grp_idx, measure = "count")
p_wgt <- compareInteractions(cellchat, show.legend = FALSE, group = grp_idx, measure = "weight")

pdf("results/cellchat_compare_overview_and_diff.pdf", width = 11, height = 7)
print(p_cnt + p_wgt)
par(mfrow = c(1, 2), xpd = TRUE)
# 用 invisible() 抑制 recordedplot 自动重放，确保两张差异网络并排占满整页
invisible(netVisual_diffInteraction(cellchat, weight.scale = TRUE))
invisible(netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight"))
par(mfrow = c(1, 1))
draw(netVisual_heatmap(cellchat) + netVisual_heatmap(cellchat, measure = "weight"))
# 统一最大边权尺度后再并排绘制各条件圈图，确保线宽可跨组直读
w_max <- getMaxWeight(cc_list, attribute = c("idents", "count"))
par(mfrow = c(1, length(cc_list)), xpd = TRUE)
for (nm in names(cc_list)) {
  invisible(netVisual_circle(
    cc_list[[nm]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = w_max[2],
    edge.width.max = 12,
    title.name = paste("Count -", nm)
  ))
}
dev.off()

p_rank_stack <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
p_rank_split <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
ggsave("results/cellchat_compare_ranknet.pdf", plot = p_rank_stack + p_rank_split, width = 10, height = 6, units = "in")
