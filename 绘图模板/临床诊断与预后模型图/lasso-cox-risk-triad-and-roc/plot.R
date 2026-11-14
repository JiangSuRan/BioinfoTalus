## LASSO-Cox 预后特征系数、有序风险评分三联热图与多时点独立验证 ROC 面板
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)

res <- readRDS("results/lasso_cox_model_results.rds")

sig_df <- res$signature_coef %>%
  mutate(
    Direction = ifelse(Coef >= 0, "Risk", "Protective"),
    Gene = reorder(Gene, Coef)
  )
gene_order <- rev(levels(sig_df$Gene))

# 面板 A：LASSO 非零回归权重分布
p_coef <- ggplot(sig_df, aes(x = Coef, y = Gene, fill = Direction)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  scale_fill_manual(values = c("Risk" = "#C0392B", "Protective" = "#2B5C8F")) +
  labs(title = "LASSO-Cox Signature Coefficients", x = "Coefficient", y = NULL, fill = NULL) +
  theme_classic(base_size = 10.5) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

# 面板 D：独立测试集 1/3/5 年时间依赖 ROC 曲线与内嵌 AUC 图例
roc_labels <- res$test_roc %>%
  distinct(Timepoint, AUC) %>%
  mutate(Label = sprintf("%s (AUC = %.3f)", Timepoint, AUC))
roc_df <- res$test_roc %>%
  left_join(roc_labels, by = "Timepoint")

p_roc <- ggplot(roc_df, aes(x = FP, y = TP, color = Label)) +
  geom_line(linewidth = 0.95) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Time-Dependent ROC (Test Cohort)",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.68, 0.22),
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(face = "bold", size = 11)
  )

# 按验证集风险评分升序排列病人轴，联动绘制三联分布图与热图
test_ord <- res$test_scores %>%
  arrange(RiskScore) %>%
  mutate(
    PatientIndex = row_number(),
    StatusLabel = factor(ifelse(OS == 1, "Dead", "Alive"), levels = c("Alive", "Dead"))
  )

p_risk <- ggplot(test_ord, aes(x = PatientIndex, y = RiskScore, color = RiskGroup)) +
  geom_point(size = 1.6, alpha = 0.85) +
  geom_hline(yintercept = res$risk_cutoff, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Low" = "#2B5C8F", "High" = "#C0392B")) +
  labs(title = "Test Cohort Risk Stratification", x = NULL, y = "Risk Score", color = "Risk") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 11)
  )

p_status <- ggplot(test_ord, aes(x = PatientIndex, y = OS.time / 365, color = StatusLabel)) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_color_manual(values = c("Alive" = "#2B5C8F", "Dead" = "#C0392B")) +
  labs(x = NULL, y = "OS (Years)", color = "Status") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")

# 行基因按与面板 A 相同的权重降序排列并截断 Z-score，使热图上下分区直觉一致
expr_ord <- t(res$test_expr_sig[test_ord$Sample, gene_order, drop = FALSE])
expr_z <- t(scale(t(expr_ord)))
expr_z[expr_z > 2.5] <- 2.5
expr_z[expr_z < -2.5] <- -2.5

ha_top <- HeatmapAnnotation(
  Risk = test_ord$RiskGroup,
  col = list(Risk = c("Low" = "#2B5C8F", "High" = "#C0392B")),
  show_annotation_name = FALSE
)
ht <- Heatmap(
  expr_z,
  name = "Z-score",
  col = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B")),
  heatmap_legend_param = list(at = c(-2.5, 0, 2.5)),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9),
  top_annotation = ha_top
)
p_heat <- wrap_elements(full = grid.grabExpr(draw(ht, merge_legend = TRUE)))

p_left <- p_coef / p_roc + plot_layout(heights = c(0.95, 1.15))
p_right <- p_risk / p_status / p_heat + plot_layout(heights = c(0.9, 0.9, 1.5))
p_all <- (p_left | p_right) + plot_layout(widths = c(1, 1.45))

ggsave("results/lasso_cox_risk_triad_roc.pdf", p_all, width = 13, height = 8.5, units = "in")
