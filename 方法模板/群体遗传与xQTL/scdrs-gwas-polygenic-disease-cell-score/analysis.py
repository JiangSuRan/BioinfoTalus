## scDRS 将 GWAS 多基因风险基因集映射到单细胞尺度估算疾病关联得分与异质性检验
import pandas as pd
import scanpy as sc
import scdrs

adata = sc.read_h5ad("Data/scrna_preprocessed.h5ad")
df_cov = pd.read_csv("Data/cell_covariates.tsv", sep="\t", index_col=0)
dict_gs = scdrs.util.load_gs("Data/trait_magma.gs")

trait_name = list(dict_gs.keys())[0]
gene_list, gene_weight = dict_gs[trait_name]

# 回归技术协变量后按基因均值-方差分层匹配 1000 组对照基因集，估算单细胞疾病关联 Z 分数与 P 值
scdrs.preprocess(adata, cov=df_cov)
df_score = scdrs.score_cell(
    data=adata,
    gene_list=gene_list,
    gene_weight=gene_weight,
    ctrl_match_key="mean_var",
    n_ctrl=1000,
    weight_opt="vs",
    return_ctrl_raw_score=False,
    return_ctrl_norm_score=True,
    random_seed=42,
)

# 结合细胞近邻图检验目标亚群内部的疾病相关性整体富集度与群内异质性（MC 置换）
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
group_res = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=df_score,
    group_cols=["cell_type"],
)

adata.obs["scdrs_norm_score"] = df_score["norm_score"]
adata.obs["scdrs_pval"] = df_score["pval"]
df_score.to_csv("results/scdrs_cell_scores.csv")
group_res["cell_type"].to_csv("results/scdrs_group_association_heterogeneity.csv")
adata.write_h5ad("results/scdrs_scored.h5ad", compression="gzip")
