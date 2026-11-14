## 本地 eQTL 汇总统计 PLINK 去连锁、弱工具剔除、结局 VCF 检索与两样本 MR 敏感性分析
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MRPRESSO)

eqtl_dt <- fread("Data/eqtlgen_target_subset.tsv")
af_dt <- fread("Data/allele_frequency.tsv")
outcome_dt <- fread("Data/outcome_snp_lookup.tsv")
plink_bin <- "Data/plink"
ld_ref_prefix <- "Data/1000G_EUR_phase3"

# 从 Zscore 与样本量 N 推导效应量 beta 与标准误 se，计算 F 统计量过滤弱工具变量（F > 10）
sig_eqtl <- eqtl_dt %>%
  filter(!is.na(Pvalue), Pvalue < 5e-8) %>%
  left_join(af_dt %>% select(SNP, MAF), by = "SNP") %>%
  mutate(
    maf_use = ifelse(is.na(MAF) | MAF <= 0 | MAF >= 1, 0.25, MAF),
    se = 1 / sqrt(2 * maf_use * (1 - maf_use) * (NrSamples + Zscore^2)),
    beta = Zscore * se,
    R2 = (2 * beta^2 * maf_use * (1 - maf_use)) / (2 * beta^2 * maf_use * (1 - maf_use) + se^2 * 2 * NrSamples * maf_use * (1 - maf_use)),
    F_stat = (R2 * (NrSamples - 2)) / (1 - R2)
  ) %>%
  filter(F_stat > 10)

# 按基因调用本地 PLINK 严格去除连锁不平衡（r2 < 0.001, 100 kb）保留独立工具变量
clump_list <- vector("list", length(unique(sig_eqtl$gene)))
names(clump_list) <- unique(sig_eqtl$gene)
for (g in names(clump_list)) {
  sub_dt <- sig_eqtl %>% filter(gene == g) %>% select(SNP, P = Pvalue) %>% distinct(SNP, .keep_all = TRUE)
  in_tsv <- file.path("results", paste0(g, "_clump_in.tsv"))
  out_prefix <- file.path("results", paste0(g, "_clump"))
  fwrite(sub_dt, in_tsv, sep = "\t")
  system2(plink_bin, args = c("--bfile", ld_ref_prefix, "--clump", in_tsv, "--clump-snp-field", "SNP", "--clump-field", "P", "--clump-p1", "5e-8", "--clump-r2", "0.001", "--clump-kb", "100", "--out", out_prefix))
  lead_snps <- fread(paste0(out_prefix, ".clumped"))$SNP
  clump_list[[g]] <- sig_eqtl %>% filter(gene == g, SNP %in% lead_snps)
}
iv_dt <- bind_rows(clump_list)

exp_fmt <- iv_dt %>%
  transmute(
    SNP = SNP,
    id.exposure = gene,
    exposure = gene,
    beta.exposure = beta,
    se.exposure = se,
    pval.exposure = Pvalue,
    effect_allele.exposure = AssessedAllele,
    other_allele.exposure = OtherAllele,
    eaf.exposure = maf_use,
    samplesize.exposure = NrSamples
  )

out_fmt <- outcome_dt %>%
  transmute(
    SNP = SNP,
    id.outcome = "Outcome",
    outcome = "Outcome",
    beta.outcome = beta,
    se.outcome = se,
    pval.outcome = pval,
    effect_allele.outcome = effect_allele,
    other_allele.outcome = other_allele,
    eaf.outcome = eaf
  )

harm <- harmonise_data(exp_fmt, out_fmt, action = 3) %>%
  filter(mr_keep, is.na(pval.outcome) | pval.outcome >= 5e-8)

het_res <- mr_heterogeneity(harm)
# 根据 Cochran Q 异质性显著性与否自适应选择固定效应或随机效应 IVW 主估计器
ivw_method <- ifelse(any(het_res$Q_pval < 0.05, na.rm = TRUE), "mr_ivw_mre", "mr_ivw_fe")
mr_res <- mr(harm, method_list = c(ivw_method, "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")) %>%
  generate_odds_ratios()
pleio_res <- mr_pleiotropy_test(harm)
presso_res <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  data = as.data.frame(harm),
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  SignifThreshold = 0.05,
  NbDistribution = 1000,
  seed = 123
)
loo_res <- mr_leaveoneout(harm)
single_res <- mr_singlesnp(harm)

fwrite(harm, "results/mr_analysis_ready.tsv", sep = "\t")
fwrite(mr_res, "results/mr_results_all_methods.tsv", sep = "\t")
fwrite(het_res, "results/mr_heterogeneity.tsv", sep = "\t")
fwrite(pleio_res, "results/mr_pleiotropy.tsv", sep = "\t")
saveRDS(presso_res, "results/mr_presso_results.rds")
fwrite(loo_res, "results/mr_leaveoneout.tsv", sep = "\t")
fwrite(single_res, "results/mr_singlesnp.tsv", sep = "\t")
