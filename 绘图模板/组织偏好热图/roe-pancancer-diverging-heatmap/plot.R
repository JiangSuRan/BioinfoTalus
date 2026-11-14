## 泛癌队列焦点组织 Ro/e 发散色阶热图
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

meta <- readRDS("Data/pancancer_cell_tissue_metadata.rds")

roe_list <- list()
for (cohort_i in sort(unique(meta$cohort))) {
  sub_df <- meta %>% filter(cohort == cohort_i)
  if (!("Target" %in% sub_df$tissue_group) || sum(sub_df$tissue_group != "Target") == 0) {
    next
  }
  obs <- table(as.character(sub_df$cell_type), as.character(sub_df$tissue_group))
  exp_mat <- outer(rowSums(obs), colSums(obs)) / sum(obs)
  roe_mat <- unclass(obs) / exp_mat
  cohort_label <- paste0(cohort_i, " (N = ", length(unique(sub_df$sample_id)), ")")
  roe_list[[cohort_label]] <- data.frame(
    cohort = cohort_label,
    cell_type = rownames(roe_mat),
    roe = roe_mat[, "Target"],
    stringsAsFactors = FALSE
  )
}

roe_wide <- bind_rows(roe_list) %>%
  pivot_wider(names_from = cell_type, values_from = roe) %>%
  as.data.frame()
rownames(roe_wide) <- roe_wide$cohort
roe_wide$cohort <- NULL
roe_mat_all <- as.matrix(roe_wide)

col_fun <- colorRamp2(c(0, 0.5, 1, 1.5, 2), c("#0081a7", "#00afb9", "#fdfcdc", "#fed9b7", "#f07167"))

ht <- Heatmap(
  roe_mat_all,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  column_names_side = "bottom",
  name = "Ro/e",
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    title_position = "topcenter",
    legend_width = unit(1.0, "inch"),
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(roe_mat_all[i, j])) {
      grid.text(sprintf("%.1f", roe_mat_all[i, j]), x, y, gp = gpar(fontsize = 6))
    }
  },
  width = ncol(roe_mat_all) * unit(0.45, "inch"),
  height = nrow(roe_mat_all) * unit(0.25, "inch"),
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7)
)

pdf("results/roe_pancancer_diverging_heatmap.pdf", width = 5.5, height = 6)
draw(ht, heatmap_legend_side = "bottom")
dev.off()
