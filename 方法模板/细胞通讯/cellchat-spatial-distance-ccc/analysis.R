## CellChat v2 空间距离约束细胞通讯推断与空间信号分布可视化
library(Seurat)
library(CellChat)
library(dplyr)

visium_obj <- readRDS("Data/visium_annotated.rds")

data_input <- LayerData(visium_obj, assay = "Spatial", layer = "data")
meta_df <- data.frame(
  labels = Idents(visium_obj),
  samples = factor(visium_obj$orig.ident),
  row.names = colnames(visium_obj)
)
spatial_locs <- GetTissueCoordinates(visium_obj, scale = NULL)
coord_mat <- as.matrix(spatial_locs[, c("x", "y")])

# 换算像素到微米尺度比值与 Visium 点径参数，以物理距离截断非邻近细胞间的旁分泌通讯
conversion_factor <- 65 / 120
spatial_factors <- data.frame(ratio = conversion_factor, tol = 65 / 2)

cellchat <- createCellChat(
  object = data_input,
  meta = meta_df,
  group.by = "labels",
  datatype = "spatial",
  coordinates = coord_mat,
  spatial.factors = spatial_factors
)

cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.1,
  distance.use = TRUE,
  interaction.range = 250,
  scale.distance = 0.01,
  contact.dependent = TRUE,
  contact.range = 100
)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "results/cellchat_spatial_ccc.rds")
write.csv(subsetCommunication(cellchat), "results/cellchat_spatial_lr_pairs.csv", row.names = FALSE)
