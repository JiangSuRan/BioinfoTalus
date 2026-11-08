## CellChat 多条件对象合并、通讯数量/强度比较与信号通路 rankNet
library(CellChat)
library(patchwork)

cellchat.NL <- readRDS("Data/cellchat_WT12_final.RDS")
cellchat.LS <- readRDS("Data/cellchat_DARE12_final.RDS")
cellchat.NL <- updateCellChat(cellchat.NL)
cellchat.LS <- updateCellChat(cellchat.LS)
object.list <- list(WT12 = cellchat.NL, DARE12 = cellchat.LS)

# 合并多条件对象（对齐细胞类型与信号基因）
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 通讯数量与强度组间条形图
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")

# 信号通路信息流比较（堆叠与威尔科克森检验版）
gg3 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg4 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)

pdf("results/cellchat_comparison.pdf", width = 12, height = 8)
print(gg1 + gg2)
print(gg3 + gg4)
dev.off()

saveRDS(cellchat, "results/cellchatCombined_conditions.RDS")
