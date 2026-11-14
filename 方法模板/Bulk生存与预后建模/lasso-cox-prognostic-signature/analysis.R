## 分层划分队列的 LASSO-Cox 预后风险签名构建与独立验证时间依赖 ROC 评估
library(survival)
library(glmnet)
library(caret)
library(survivalROC)

cohort <- readRDS("Data/cohort_survival_expr.rds")
cohort <- cohort[!is.na(cohort$OS.time) & !is.na(cohort$OS) & cohort$OS.time > 0, , drop = FALSE]
gene_cols <- setdiff(colnames(cohort), c("OS.time", "OS"))

# 按删失结局状态分层划分 70% 训练与 30% 验证队列，防止结局分布漂移与测试泄露
set.seed(123)
train_idx <- createDataPartition(cohort$OS, p = 0.7, list = FALSE)[, 1]
train_df <- cohort[train_idx, , drop = FALSE]
test_df <- cohort[-train_idx, , drop = FALSE]

# 仅在训练集执行单因素 Cox 预筛选，避免测试集信息泄露到特征集合
uni_p <- sapply(gene_cols, function(g) {
  fit <- coxph(Surv(OS.time, OS) ~ train_df[[g]], data = train_df)
  summary(fit)$coefficients[1, "Pr(>|z|)"]
})
cand_genes <- names(uni_p)[!is.na(uni_p) & uni_p < 0.05]

x_train <- as.matrix(train_df[, cand_genes, drop = FALSE])
y_train <- Surv(train_df$OS.time, train_df$OS)

# 10 折交叉验证严格使用 alpha = 1 的 LASSO 惩罚，并在 lambda.min 提取非零签名系数
set.seed(123)
cv_fit <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 1, nfolds = 10, cox.ties = "efron")
coef_vec <- as.numeric(coef(cv_fit, s = "lambda.min"))
names(coef_vec) <- colnames(x_train)
active_coef <- coef_vec[coef_vec != 0]
sig_genes <- names(active_coef)

# 风险评分阈值严格取自训练队列中位数并固定迁移到独立测试集
train_risk <- as.numeric(x_train[, sig_genes, drop = FALSE] %*% active_coef)
risk_cutoff <- median(train_risk)

x_test <- as.matrix(test_df[, sig_genes, drop = FALSE])
test_risk <- as.numeric(x_test %*% active_coef)
test_scores <- data.frame(
  Sample = rownames(test_df),
  OS.time = test_df$OS.time,
  OS = test_df$OS,
  RiskScore = test_risk,
  RiskGroup = factor(ifelse(test_risk > risk_cutoff, "High", "Low"), levels = c("Low", "High")),
  stringsAsFactors = FALSE
)

# 独立验证集时间依赖 ROC（1/3/5 年，NNE 平滑保证敏感度单调有界于 [0, 1]）
eval_times <- c("1-Year" = 365, "3-Year" = 1095, "5-Year" = 1825)
nne_span <- 0.25 * nrow(test_scores)^(-0.20)
roc_list <- lapply(names(eval_times), function(lbl) {
  r <- survivalROC(
    Stime = test_scores$OS.time,
    status = test_scores$OS,
    marker = test_scores$RiskScore,
    predict.time = eval_times[[lbl]],
    method = "NNE",
    span = nne_span
  )
  data.frame(Timepoint = lbl, FP = r$FP, TP = r$TP, AUC = r$AUC, stringsAsFactors = FALSE)
})
test_roc <- do.call(rbind, roc_list)

km_diff <- survdiff(Surv(OS.time, OS) ~ RiskGroup, data = test_scores)
test_metrics <- data.frame(
  Timepoint = names(eval_times),
  AUC = sapply(roc_list, function(d) d$AUC[1]),
  LogRank_P = pchisq(km_diff$chisq, df = 1, lower.tail = FALSE),
  row.names = NULL
)
sig_df <- data.frame(Gene = sig_genes, Coef = unname(active_coef), stringsAsFactors = FALSE)

model_bundle <- list(
  cv_fit = cv_fit,
  risk_cutoff = risk_cutoff,
  signature_coef = sig_df,
  test_scores = test_scores,
  test_expr_sig = x_test,
  test_roc = test_roc,
  test_metrics = test_metrics
)
saveRDS(model_bundle, "results/lasso_cox_model_results.rds")
write.csv(sig_df, "results/lasso_cox_signature_genes.csv", row.names = FALSE)
write.csv(test_metrics, "results/lasso_cox_test_metrics.csv", row.names = FALSE)
