## scRepertoire 免疫组库克隆动态、重叠热图与亚群共享弦图
library(scRepertoire)
library(Seurat)
library(ggplot2)
library(circlize)
library(scales)
library(patchwork)

combined_tcr <- readRDS("Data/combined_tcr.rds")
sc_vdj <- readRDS("Data/seurat_vdj.rds")

p_quant <- clonalQuant(combined_tcr, cloneCall = "strict", scale = TRUE) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75))

p_homeo <- clonalHomeostasis(combined_tcr, cloneCall = "gene")
p_ovl <- clonalOverlap(combined_tcr, cloneCall = "strict", method = "morisita")

p_overlay <- clonalOverlay(
  sc_vdj,
  reduction = "umap",
  cut.category = "clonalFrequency",
  cutpoint = 1,
  bins = 10,
  facet.by = "orig.ident"
) + guides(color = "none")

p_occ <- clonalOccupy(sc_vdj, x.axis = "seurat_clusters", proportion = TRUE, label = FALSE)
p_alluv <- alluvialClones(
  sc_vdj,
  cloneCall = "gene",
  y.axes = c("orig.ident", "seurat_clusters"),
  color = "orig.ident"
)

p_top <- (p_quant | p_homeo) / (p_ovl | p_occ)
ggsave("results/screpertoire_clonal_suite.pdf", plot = (p_top / (p_overlay | p_alluv)), width = 12, height = 11, units = "in")

circles <- getCirclize(sc_vdj, cloneCall = "gene", group.by = "seurat_clusters", proportion = TRUE)
cluster_lvls <- levels(factor(sc_vdj$seurat_clusters))
grid_cols <- setNames(hue_pal()(length(cluster_lvls)), cluster_lvls)

pdf("results/screpertoire_cluster_chord.pdf", width = 7, height = 7)
chordDiagram(circles, self.link = 1, grid.col = grid_cols, directional = 1, direction.type = "arrows", link.arr.type = "big.arrow")
dev.off()
