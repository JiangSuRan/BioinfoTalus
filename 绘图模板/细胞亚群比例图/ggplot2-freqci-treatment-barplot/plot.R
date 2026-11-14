## 多项分布 95% 置信区间与卡方检验亚群治疗反应比例柱状图
library(REdaS)
library(dplyr)
library(forcats)
library(ggplot2)

meta <- readRDS("Data/cell_treatment_metadata.rds")
meta$group <- as_factor(meta$group)

ci_list <- list()
for (grp in levels(meta$group)) {
  sub_vec <- factor(meta$cell_type[meta$group == grp], levels = sort(unique(meta$cell_type)))
  ci_mat <- as.data.frame.array(print(freqCI(sub_vec, level = 0.95), percent = TRUE, digits = 3))
  colnames(ci_mat)[c(1, 3)] <- c("lower", "upper")
  ci_mat$cell_type <- rownames(ci_mat)
  ci_mat$group <- grp
  ci_list[[grp]] <- ci_mat
}

tab_ci <- bind_rows(ci_list)
tab_ci$group <- factor(tab_ci$group, levels = levels(meta$group))

count_mat <- as.matrix(table(meta$group, meta$cell_type))
chi_p <- apply(count_mat, 2, function(x) chisq.test(x)$p.value)
chi_df <- data.frame(
  cell_type = colnames(count_mat),
  p_value = chi_p,
  p_bonferroni = pmin(chi_p * ncol(count_mat), 1)
)
write.csv(chi_df, "results/celltype_proportion_chisq.csv", row.names = FALSE)

p <- ggplot(tab_ci, aes(x = group, y = Estimate, fill = cell_type)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.05)) +
  facet_grid(cols = vars(cell_type)) +
  labs(x = NULL, y = "Percentage of cells (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )

ggsave(filename = "results/treatment_freqci_proportion_barplot.pdf", plot = p, width = 11, height = 6, units = "in")
