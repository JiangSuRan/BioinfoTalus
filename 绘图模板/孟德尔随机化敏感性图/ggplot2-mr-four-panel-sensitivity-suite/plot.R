## 两样本孟德尔随机化四联敏感性诊断图（散点/单SNP森林/漏斗/留一法）
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

dat <- fread("Data/mr_analysis_ready.tsv") %>% as.data.frame()
res_all <- fread("Data/mr_results_all_methods.tsv") %>% as.data.frame()
res_single <- fread("Data/mr_singlesnp.tsv") %>% as.data.frame()
res_loo <- fread("Data/mr_leaveoneout.tsv") %>% as.data.frame()

# 当暴露效应方向与正向等位基因编码相反时翻转符号至第一象限，保持散点与斜率直观一致
sign_flip <- ifelse(dat$beta.exposure < 0, -1, 1)
dat <- dat %>%
  mutate(
    beta.exposure.pos = abs(beta.exposure),
    beta.outcome.pos = beta.outcome * sign_flip
  )

method_pal <- c(
  "Inverse variance weighted" = "#1F78B4",
  "MR Egger" = "#E31A1C",
  "Weighted median" = "#33A02C",
  "Simple mode" = "#FF7F00",
  "Weighted mode" = "#6A3D9A"
)

slope_df <- res_all %>%
  filter(method %in% names(method_pal)) %>%
  mutate(intercept = if_else(method == "MR Egger" & "egger_intercept" %in% names(res_all), egger_intercept, 0))

p_scatter <- ggplot(dat, aes(x = beta.exposure.pos, y = beta.outcome.pos)) +
  geom_errorbar(aes(ymin = beta.outcome.pos - 1.96 * se.outcome, ymax = beta.outcome.pos + 1.96 * se.outcome), width = 0, colour = "grey65") +
  geom_errorbarh(aes(xmin = beta.exposure.pos - 1.96 * se.exposure, xmax = beta.exposure.pos + 1.96 * se.exposure), height = 0, colour = "grey65") +
  geom_point(size = 2.2, colour = "#2C3E50") +
  geom_abline(data = slope_df, aes(slope = b, intercept = intercept, colour = method), linewidth = 0.85) +
  scale_colour_manual(values = method_pal) +
  labs(title = "A. MR Method Slopes", x = "SNP effect on exposure", y = "SNP effect on outcome", colour = "MR Method") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

single_df <- res_single %>%
  mutate(
    ci_low = b - 1.96 * se,
    ci_high = b + 1.96 * se,
    is_summary = grepl("^All -", SNP),
    SNP = factor(SNP, levels = rev(SNP))
  )

p_forest <- ggplot(single_df, aes(x = b, y = SNP, colour = is_summary)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
  geom_point(aes(size = is_summary)) +
  scale_colour_manual(values = c("FALSE" = "#34495E", "TRUE" = "#E64B35"), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 2.8), guide = "none") +
  labs(title = "B. Single-SNP Wald Ratio", x = "Causal effect (beta +/- 95% CI)", y = "") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

funnel_df <- single_df %>% filter(!is_summary) %>% mutate(inv_se = 1 / se)
ivw_summary <- res_all %>% filter(method %in% c("Inverse variance weighted", "MR Egger"))

p_funnel <- ggplot(funnel_df, aes(x = b, y = inv_se)) +
  geom_point(size = 2.1, colour = "#2C3E50") +
  geom_vline(data = ivw_summary, aes(xintercept = b, colour = method), linewidth = 0.8) +
  scale_colour_manual(values = method_pal[c("Inverse variance weighted", "MR Egger")]) +
  labs(title = "C. Funnel Symmetry", x = expression(beta[IV]), y = expression(1 / SE[IV]), colour = "Estimator") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

loo_df <- res_loo %>%
  mutate(
    ci_low = b - 1.96 * se,
    ci_high = b + 1.96 * se,
    is_all = SNP == "All",
    SNP = factor(SNP, levels = rev(SNP))
  )

p_loo <- ggplot(loo_df, aes(x = b, y = SNP, colour = is_all)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
  geom_point(aes(size = is_all)) +
  scale_colour_manual(values = c("FALSE" = "#3B4992", "TRUE" = "#E64B35"), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 2.8), guide = "none") +
  labs(title = "D. Leave-One-Out Sensitivity", x = "IVW beta after excluding SNP", y = "") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())

p_suite <- (p_scatter | p_forest) / (p_funnel | p_loo)

ggsave(filename = "results/mr_four_panel_sensitivity.pdf", plot = p_suite, width = 12.5, height = 10, units = "in")
