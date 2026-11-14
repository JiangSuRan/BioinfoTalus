## NicheNet 优先配体表达分布与配体-靶基因/受体调控热图
library(Seurat)
library(nichenetr)
library(ggplot2)
library(patchwork)

sc_obj <- readRDS("Data/seurat_obj.rds")
nn_out <- readRDS("Data/nichenet_output.rds")
Idents(sc_obj) <- "celltype"

p_lig <- DotPlot(sc_obj, features = rev(nn_out$top_ligands), cols = "RdYlBu") +
  RotatedAxis() +
  ggtitle("Prioritized Ligands across Cell Types")

p_lt <- nn_out$ligand_target_heatmap +
  scale_fill_gradient2(low = "whitesmoke", high = "royalblue") +
  labs(title = "Ligand-Target Regulatory Potential")

p_lr <- nn_out$ligand_receptor_heatmap +
  labs(title = "Ligand-Receptor Interaction")

p_all <- p_lig / (p_lt | p_lr)
ggsave("results/nichenet_ligand_target_receptor.pdf", plot = p_all, width = 12, height = 10, units = "in")
