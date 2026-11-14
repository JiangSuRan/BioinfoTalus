## SPATA2 组织边界外点剔除、SPARK-X 空间基因筛选与定量空间轨迹梯度推断
library(SPATA2)
library(dplyr)

spata_obj <- readRDS("Data/spata2_object.rds")

# 识别组织形态学轮廓并剔除游离空间离群 Spot，避免边界假性梯度干扰
spata_obj <- identifyTissueOutline(spata_obj)
spata_obj <- removeSpatialOutliers(spata_obj, verbose = FALSE)
spata_obj <- runSPARKX(spata_obj)

# 定义方向性空间坐标轨迹轴并筛选沿该组织梯度呈现单调或线性变化模式基因（STS）
spata_obj <- addSpatialTrajectory(
  object = spata_obj,
  id = "Hypoxia_Axis",
  start = c(x = 150, y = 200),
  end = c(x = 650, y = 550),
  width = 120
)

sts_res <- spatialTrajectoryScreening(
  object = spata_obj,
  id = "Hypoxia_Axis",
  variables = getSparkxGenes(spata_obj, threshold_pval = 0.01)
)

saveRDS(spata_obj, "results/spata2_trajectory_screened.rds")
write.csv(sts_res@results, "results/spata2_trajectory_genes.csv", row.names = FALSE)
