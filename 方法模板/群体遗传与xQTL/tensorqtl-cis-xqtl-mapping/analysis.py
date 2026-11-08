## tensorQTL 单细胞 Pseudobulk 顺式表达/染色质可及性数量性状位点（cis-eQTL / cis-caQTL）GPU 映射
import pandas as pd
import torch
import genotypeio
import tensorqtl
from tensorqtl import cis

# 逆正态变换表型（bed.gz）与协变量按样本对齐
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed("Data/expression_bed.bed.gz")
covariates_df = pd.read_csv("Data/covariates.csv", index_col=0)
phenotype_df = phenotype_df[covariates_df.index]

# PEER 因子数按样本量自适应分档，避免小样本过拟合
sample_size = covariates_df.shape[0]
if sample_size <= 150:
    n_peer = 15
elif sample_size <= 250:
    n_peer = 30
else:
    n_peer = 35
peer_cols = ["pf" + str(i) for i in range(1, n_peer + 1)]
covariates_df = covariates_df.loc[:, ["Age", "sex", "PC1", "PC2"] + peer_cols]

# PLINK 剂量基因型加载，MAF 过滤在 map_cis 内执行
reader = genotypeio.PlinkReader("Data/genotype/CIMA")
genotype_df = reader.load_genotypes()
variant_df = reader.variant_df

# beta 逼近置换检验：每表型取最优 variant 的经验 p 值，再做全局 q 值校正
cis_df = cis.map_cis(
    genotype_df, variant_df, phenotype_df, phenotype_pos_df,
    covariates_df=covariates_df, maf_threshold=0.05, seed=2025,
)
cis_df = tensorqtl.calculate_qvalues(cis_df)
cis_df.to_csv("results/cis_qtl_permutation.csv")

# 名义关联（每表型全部 variant）以 parquet 分文件输出
cis.map_nominal(
    genotype_df, variant_df, phenotype_df, phenotype_pos_df, "all",
    covariates_df=covariates_df, maf_threshold=0.05,
    output_dir="results/",
)
