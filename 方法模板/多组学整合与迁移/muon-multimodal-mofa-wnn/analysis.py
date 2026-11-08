## muon 多视图 LIANA 通讯视图构建与 MOFA+ 因子分解
import gc
import os

import liana as li
import muon as mu
import scanpy as sc

adata = sc.read_h5ad("Data/ccc_annotated.h5ad")

score_key = "liana_res"
sample_key = "sample"
condition_key = "condition"
n_factors = 10

# LIANA 逐样本通讯得分 → views（样本 × LR 矩阵矩阵化）
mdata = li.multi.lrs_to_views(
    adata,
    sample_key=sample_key,
    score_key=score_key,
    inverse_fun=None,
    obs_keys=[condition_key],
    lr_prop=0.33,
    lrs_per_sample=3,
    lrs_per_view=15,
    samples_per_view=5,
    min_variance=0,
    lr_fill=0,
    verbose=True,
    uns_key=score_key,
).copy()

os.makedirs("results/mofa", exist_ok=True)

# MOFA+ 多视图因子分解
mu.tl.mofa(
    mdata,
    use_obs="union",
    convergence_mode="medium",
    n_factors=n_factors,
    outfile="results/mofa/liana_mofa.hdf5",
    seed=1337,
)

# 因子得分提取并与条件关联
factor_scores = li.ut.get_factor_scores(mdata, obsm_key="X_mofa").copy()
factor_scores.to_csv("results/mofa/factor_scores.csv")

gc.collect()
