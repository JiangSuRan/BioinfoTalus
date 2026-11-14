## NicheNet 发送-接收亚群配体优先级与靶基因调控推断
library(nichenetr)
library(Seurat)
library(dplyr)

sc_obj <- readRDS("Data/seurat_obj.rds")
ligand_target_matrix <- readRDS("Data/ligand_target_matrix.rds")
lr_network <- readRDS("Data/lr_network.rds")
weighted_networks <- readRDS("Data/weighted_networks.rds")

Idents(sc_obj) <- "celltype"
cond_levels <- levels(factor(sc_obj$condition))

# 将发送群表达配体与接收群条件间差异基因联合匹配先验调控网络，排序因果解释力最高的配体
nichenet_out <- nichenet_seuratobj_aggregate(
  seurat_obj = sc_obj,
  receiver = levels(Idents(sc_obj))[1],
  sender = levels(Idents(sc_obj))[-1],
  condition_colname = "condition",
  condition_oi = cond_levels[2],
  condition_reference = cond_levels[1],
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)

saveRDS(nichenet_out, "results/nichenet_output.rds")
write.table(nichenet_out$ligand_activities, "results/nichenet_ligand_activities.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
