## ggsurvfit 临床级生存曲线（survfit2 + 置信区间 + 风险表）
library(survival)
library(ggsurvfit)
library(ggprism)
library(dplyr)

surv_df <- read.csv("Data/survival_signature.csv", stringsAsFactors = FALSE)

# 按签名中位 cutoff 二分高低组
med_val <- median(surv_df$signature_score, na.rm = TRUE)
surv_df <- surv_df %>%
  mutate(median_cutoff = case_when(
    signature_score < med_val ~ "Low",
    signature_score >= med_val ~ "High"
  ))

# survfit2 + ggsurvfit 链：置信区间 + p 值 + 风险表 + 删失标记 + prism 主题
p <- survfit2(Surv(pfscdy, censpfs) ~ median_cutoff, data = surv_df) %>%
  ggsurvfit(size = 2) +
  add_confidence_interval() +
  add_pvalue() +
  add_risktable() +
  labs(x = "Days", y = "Progression Free Survival") +
  add_censor_mark(size = 3) +
  theme_prism() +
  scale_color_manual(values = c("Low" = "#4275AE", "High" = "#D69F4C"))

ggsave(filename = "results/ggsurvfit_km_curve.pdf", plot = p, width = 6, height = 6, units = "in")
