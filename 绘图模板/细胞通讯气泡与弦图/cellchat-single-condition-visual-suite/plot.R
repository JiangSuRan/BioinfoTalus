## 单条件 CellChat 通讯网络、配体受体气泡及信号发送接收角色全景图
library(CellChat)
library(patchwork)
library(ComplexHeatmap)

cellchat <- readRDS("Data/cellchat_obj.rds")
group_size <- as.numeric(table(cellchat@idents))
target_pathway <- cellchat@netP$pathways[1]
n_groups <- nlevels(cellchat@idents)

pdf("results/cellchat_single_network_overview.pdf", width = 11, height = 8)
par(mfrow = c(1, 2), xpd = TRUE)
# 用 invisible() 抑制 recordedplot 自动重放，避免 mfrow 子图被重置为单页左上角
invisible(netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions"))
invisible(netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength"))
par(mfrow = c(1, 1), xpd = TRUE)
invisible(netVisual_aggregate(cellchat, signaling = target_pathway, layout = "chord", vertex.weight = group_size))
draw(netVisual_heatmap(cellchat, signaling = target_pathway, color.heatmap = "Reds"))
print(netAnalysis_contribution(cellchat, signaling = target_pathway))
dev.off()

p_bubble <- netVisual_bubble(cellchat, sources.use = seq_len(min(2, n_groups)), targets.use = seq_len(n_groups), remove.isolate = FALSE)
p_role <- netAnalysis_signalingRole_scatter(cellchat)
ht_out <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht_in <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

pdf("results/cellchat_single_bubble_and_roles.pdf", width = 10, height = 8)
print(p_bubble)
print(p_role)
draw(ht_out + ht_in)
dev.off()
