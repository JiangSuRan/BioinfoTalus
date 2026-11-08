## 临床队列单因素/多因素 Cox 比例风险回归（coxph）与预后风险评分建模
library(survival)

clinical_dat <- read.csv("Data/clinical_survival.csv", stringsAsFactors = FALSE)

# 生存对象：随访时间（月）与事件结局
clinical_dat$SurvObj <- with(clinical_dat, Surv(Time_to_recurrence_months, recurrence))

# 逐变量单因素 Cox：各变量独立效应的粗筛与估计
uni_vars <- c("Gender", "Tumor_size_cm", "Grade", "Stage", "TumorMacrophageGSEA")
uni_res <- do.call(rbind, lapply(uni_vars, function(v) {
  fit <- coxph(as.formula(paste("SurvObj ~", v)), data = clinical_dat)
  s <- summary(fit)
  data.frame(
    variable = v,
    HR = round(s$conf.int[1, 1], 2),
    HR.lower = round(s$conf.int[1, 3], 2),
    HR.upper = round(s$conf.int[1, 4], 2),
    coef = s$coefficients[1, 1],
    Pval = round(s$coefficients[1, 5], 4)
  )
}))
uni_res$adj.Pval <- p.adjust(uni_res$Pval, method = "BH")

# 多因素模型：候选临床变量与分子评分联合建模
multi_fit <- coxph(as.formula(paste("SurvObj ~", paste(uni_vars, collapse = " + "))), data = clinical_dat)
s <- summary(multi_fit)
multi_res <- data.frame(
  variable = rownames(s$conf.int),
  HR = round(s$conf.int[, 1], 2),
  HR.lower = round(s$conf.int[, 3], 2),
  HR.upper = round(s$conf.int[, 4], 2),
  coef = s$coefficients[, 1],
  Pval = round(s$coefficients[, 5], 4)
)
multi_res$adj.Pval <- p.adjust(multi_res$Pval, method = "BH")

# 预后风险评分：多因素模型系数对各协变量的线性组合
clinical_dat$risk_score <- predict(multi_fit, type = "lp")

cox_res <- rbind(
  cbind(uni_res, model = "univariate"),
  cbind(multi_res, model = "multivariate")
)
write.csv(cox_res, "results/cox_regression_results.csv", row.names = FALSE)
saveRDS(clinical_dat, "results/clinical_with_risk_score.rds")
