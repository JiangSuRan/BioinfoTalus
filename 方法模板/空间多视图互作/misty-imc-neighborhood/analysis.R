## MISTy 成像质谱流式（IMC）单细胞空间邻域互作推断
library(mistyR)
library(dplyr)

expr <- readRDS("Data/imc_expr.rds")
pos <- readRDS("Data/imc_positions.rds")

# juxtaview 捕获直接接触的邻居，paraview 以距离衰减加权更广邻域（zoi 排除自身与近邻重复计数）
views <- create_initial_view(expr) %>%
  add_juxtaview(pos, neighbor.thr = 12) %>%
  add_paraview(pos, l = 100, zoi = 12)

run_misty(views, "results/misty_imc_out")

misty_res <- collect_results("results/misty_imc_out")
saveRDS(misty_res, "results/misty_imc_results.rds")
