## 人类 CellChat 通讯网络推断
library(CellChat)

cellchat <- readRDS("Data/cellchat_input.rds")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "results/cellchat_network.rds")
