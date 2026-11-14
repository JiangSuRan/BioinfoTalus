## CytoTRACE2 分化潜能光谱连续 UMAP 与亚群分化排序箱线图
library(Seurat)
library(CytoTRACE2)
library(dplyr)
library(ggplot2)
library(patchwork)

sc_obj <- readRDS("Data/seurat_obj.rds")
ct2_res <- readRDS("Data/cytotrace2_result.rds")

counts_mat <- GetAssayData(sc_obj, assay = "RNA", layer = "counts")
anno_df <- data.frame(phenotype = as.character(sc_obj@meta.data$celltype), row.names = colnames(sc_obj))
native_plots <- plotData(cytotrace2_result = ct2_res, annotation = anno_df, expression_data = counts_mat, is_seurat = FALSE)

umap_df <- Embeddings(sc_obj, reduction = "umap") %>% as.data.frame()
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
shared_cells <- intersect(rownames(umap_df), rownames(ct2_res))
umap_df <- bind_cols(umap_df[shared_cells, ], ct2_res[shared_cells, , drop = FALSE])

# 沿用高分化(红)至低分化/高干性(深紫)的多色渐变刻度以强化发育方向感
p_spec <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CytoTRACE2_Relative)) +
  geom_point(size = 1.6) +
  scale_colour_gradientn(
    colours = c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2"),
    na.value = "transparent",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),
    name = "Relative
potency",
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  labs(title = "CytoTRACE 2 Potency", x = "UMAP1", y = "UMAP2") +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

p_all <- (native_plots$Phenotype_UMAP | p_spec) / native_plots$CytoTRACE2_Boxplot_byPheno
ggsave("results/cytotrace2_potency_umap_panel.pdf", plot = p_all, width = 10, height = 8.5, units = "in")
