## MISTy 10x Visium 空间转录组多视图（Intraview / Paraview）互作建模
library(Seurat)
library(mistyR)
library(dplyr)

seurat_vs <- readRDS("Data/visium_seurat.rds")

# intraview 描述 spot 自身表达，paraview 汇总半径 l 内邻近 spot 的表达微环境
expr <- as.matrix(GetAssayData(seurat_vs, assay = "SCT", layer = "data"))
geometry <- GetTissueCoordinates(seurat_vs)

views <- create_initial_view(expr) %>%
  add_paraview(geometry, l = 10)

run_misty(views, "results/misty_visium_out")

misty_res <- collect_results("results/misty_visium_out")
saveRDS(misty_res, "results/misty_results.rds")
