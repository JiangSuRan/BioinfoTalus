## lme4 (glmer) 单细胞分辨率拟时序动态 eQTL（Genotype × Pseudotime 交互效应）广义线性混合模型检验
library(glue)
library(lme4)

celltype <- "CD4_Tn_CCR7"
gene <- "FOXP3"

# 单基因计数表（E 表达、G 基因型剂量、ptime 拟时序）与个体/细胞级协变量按行合并
phenotype <- read.csv(glue("Data/dynamic/{celltype}/{gene}.csv"), row.names = 1)
cov_df <- read.csv(glue("Data/dynamic/{celltype}_cov_df.csv"), row.names = 1)
data <- cbind(phenotype, cov_df)
data$total_counts <- scale(log(data$total_counts))

# 个体随机截距吸收供体内细胞相关；nAGQ=0 伪似然加速大规模单细胞拟合
full_model <- glmer(
  formula = E ~ G + (1 | sample) + Age + sex + total_counts + pct_counts_mt + PC1 + PC2 +
    mRNA_PC1 + mRNA_PC2 + mRNA_PC3 + mRNA_PC4 + mRNA_PC5 + ptime + G * ptime,
  family = "poisson", nAGQ = 0, data = data,
  control = glmerControl(optimizer = "nloptwrap")
)
null_model <- glmer(
  formula = E ~ G + (1 | sample) + Age + sex + total_counts + pct_counts_mt + PC1 + PC2 +
    mRNA_PC1 + mRNA_PC2 + mRNA_PC3 + mRNA_PC4 + mRNA_PC5 + ptime,
  family = "poisson", nAGQ = 0, data = data,
  control = glmerControl(optimizer = "nloptwrap")
)

# 似然比检验交互项；交互 β 与 SE 取自全模型固定效应
model_lrt <- anova(null_model, full_model)
beta_int <- fixef(full_model)["G:ptime"]
se_int <- sqrt(vcov(full_model)["G:ptime", "G:ptime"])

res <- data.frame(
  celltype = celltype,
  gene = gene,
  beta_Gxptime = beta_int,
  se_Gxptime = se_int,
  lrt_chisq = model_lrt$Chisq[2],
  lrt_pvalue = model_lrt$`Pr(>Chisq)`[2]
)
write.csv(res, "results/sc_dynamic_eqtl_lrt.csv", row.names = FALSE)
