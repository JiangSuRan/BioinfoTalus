## Kaplan-Meier 生存曲线与风险表组合图（survminer / ggplot2）
library(survival)
library(survminer)
library(ggsci)

dat <- read.csv("data/survival_metadata.csv", stringsAsFactors = FALSE)

# 生存时间换算为月并按免疫类型分组
dat$OS.time <- round(dat$OS.time / 30, 2)
dat$immuneType <- factor(dat$immuneType)

fit <- survfit(Surv(OS.time, OS) ~ immuneType, data = dat)

p <- ggsurvplot(
  fit,
  data = dat,
  size = 0.4,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  surv.median.line = "none",
  palette = get_palette("jco", length(unique(dat$immuneType))),
  legend.title = "",
  legend.labs = unique(as.character(dat$immuneType)),
  ggtheme = theme_bw()
) +
  xlab("Months")

ggsave(filename = "result/km_survival.pdf", plot = print(p), width = 6, height = 6, units = "in")
