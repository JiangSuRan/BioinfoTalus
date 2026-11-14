## METAFlux 基于亚群比例约束的单细胞稳态代谢通量平衡分析
library(METAFlux)
library(Seurat)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
data("human_blood", package = "METAFlux")

avg_res <- AverageExpression(sc_obj, assays = "RNA", group.by = "celltype")
mean_exp <- as.matrix(avg_res$RNA)

# 计算每群代谢反应活性评分（MRAS），再按与矩阵列名对齐的真实细胞比例求解联合 FBA
mras_scores <- calculate_reaction_score(data = mean_exp)
cell_prop <- prop.table(table(sc_obj$celltype))[colnames(mras_scores)] %>% as.numeric()
cell_prop <- cell_prop / sum(cell_prop)

sc_flux <- compute_sc_flux(
  num_cell = ncol(mras_scores),
  fraction = cell_prop,
  fluxscore = mras_scores,
  medium = human_blood
)

saveRDS(mras_scores, "results/metaflux_scores.rds")
saveRDS(sc_flux, "results/metaflux_sc_flux.rds")
