## 多基因诊断 Logistic 列线图、校准曲线、联合 ROC 与决策曲线（DCA）
library(dplyr)
library(reshape2)
library(rms)
library(pROC)
library(ggplot2)
library(patchwork)

expr_mat <- readRDS("Data/expr_matrix.rds")
sample_meta <- readRDS("Data/sample_group.rds")
hub_genes <- readLines("Data/hub_genes.txt")

# 将目标基因表达矩阵与二分类结局对齐，构造 rms 建模所需的 datadist 环境
model_df <- t(expr_mat[intersect(hub_genes, rownames(expr_mat)), , drop = FALSE]) %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  inner_join(sample_meta, by = "Sample") %>%
  mutate(Outcome = ifelse(Group == "Treat", 1L, 0L)) %>%
  select(Outcome, all_of(intersect(hub_genes, rownames(expr_mat))))

gene_cols <- setdiff(colnames(model_df), "Outcome")
dd <- datadist(model_df)
options(datadist = "dd")

fml <- as.formula(paste("Outcome ~", paste(gene_cols, collapse = " + ")))
fit_lrm <- lrm(fml, data = model_df, x = TRUE, y = TRUE)
model_df$Nomogram_Score <- predict(fit_lrm, type = "fitted")

nom <- nomogram(fit_lrm, fun = plogis, fun.at = c(0.05, 0.5, 0.95), lp = FALSE, funlabel = "Risk of Disease")
cal <- calibrate(fit_lrm, method = "boot", B = 500)

# rms::plot.nomogram 会重置图形布局参数，将列线图与校准曲线分别输出至独立画布以避免半幅留白
pdf("results/diagnostic_nomogram.pdf", width = 8.4, height = 5.4)
plot(nom, xfrac = 0.24, cex.axis = 0.85, cex.var = 0.95)
dev.off()

pdf("results/diagnostic_calibration.pdf", width = 6.2, height = 5.8)
par(mar = c(5, 5, 2.5, 1.8))
plot(cal, xlab = "Predicted Probability", ylab = "Observed Probability", subtitles = FALSE)
dev.off()

# 计算各单基因及联合列线图模型的 ROC 与 DeLong 95% CI
roc_vars <- c(gene_cols, "Nomogram_Score")
roc_coords <- vector("list", length(roc_vars))
auc_labels <- character(length(roc_vars))
for (i in seq_along(roc_vars)) {
  v <- roc_vars[i]
  r <- roc(model_df$Outcome, model_df[[v]], quiet = TRUE, ci = TRUE)
  ci_val <- as.numeric(r$ci)
  auc_labels[i] <- sprintf("%s (AUC = %.3f, 95%% CI %.3f-%.3f)", v, ci_val[2], ci_val[1], ci_val[3])
  roc_coords[[i]] <- data.frame(
    FPR = 1 - r$specificities,
    TPR = r$sensitivities,
    Model = auc_labels[i],
    stringsAsFactors = FALSE
  ) %>%
    arrange(FPR, TPR)
}
roc_df <- bind_rows(roc_coords) %>% mutate(Model = factor(Model, levels = auc_labels))

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, colour = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.9) +
  labs(title = "A. Diagnostic ROC Comparison", x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)", colour = "") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = c(0.62, 0.24), legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

# 按 Vickers & Elkin 定义在阈值网格上计算净获益，对比联合模型、全干预与零干预策略
thresholds <- seq(0.05, 0.95, by = 0.02)
n_total <- nrow(model_df)
dca_rows <- vector("list", length(thresholds) * (length(gene_cols) + 3L))
idx <- 1L
for (pt in thresholds) {
  w <- pt / (1 - pt)
  dca_rows[[idx]] <- data.frame(Threshold = pt, NetBenefit = mean(model_df$Outcome == 1L) - mean(model_df$Outcome == 0L) * w, Strategy = "Treat All")
  idx <- idx + 1L
  dca_rows[[idx]] <- data.frame(Threshold = pt, NetBenefit = 0, Strategy = "Treat None")
  idx <- idx + 1L
  for (v in roc_vars) {
    prob <- if (v == "Nomogram_Score") model_df[[v]] else predict(lrm(as.formula(paste("Outcome ~", v)), data = model_df), type = "fitted")
    tp <- sum(prob >= pt & model_df$Outcome == 1L)
    fp <- sum(prob >= pt & model_df$Outcome == 0L)
    nb <- (tp / n_total) - (fp / n_total) * w
    dca_rows[[idx]] <- data.frame(Threshold = pt, NetBenefit = nb, Strategy = v)
    idx <- idx + 1L
  }
}
dca_df <- bind_rows(dca_rows) %>%
  mutate(NetBenefit = pmax(NetBenefit, -0.15))

p_dca <- ggplot(dca_df, aes(x = Threshold, y = NetBenefit, colour = Strategy, linetype = Strategy == "Treat All")) +
  geom_line(linewidth = 0.85) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"), guide = "none") +
  coord_cartesian(ylim = c(-0.08, max(dca_df$NetBenefit, na.rm = TRUE) * 1.08)) +
  labs(title = "B. Decision Curve Analysis (DCA)", x = "High-Risk Threshold Probability", y = "Standardized Net Benefit", colour = "Strategy") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "right")

p_eval <- p_roc + p_dca + plot_layout(widths = c(1, 1.1))
ggsave(filename = "results/diagnostic_roc_dca_panel.pdf", plot = p_eval, width = 12, height = 5.2, units = "in")
