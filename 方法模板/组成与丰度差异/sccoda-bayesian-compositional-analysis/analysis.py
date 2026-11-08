## scCODA 层次贝叶斯 Dirichlet-multinomial 组成分析（HMC 采样与参考类型循环）
import pickle
from datetime import datetime

import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod

cell_count = pd.read_csv("Data/whole_meta.csv")
split_key = "class"
split_value = ["excitatory", "inhibitory"]
random_effect = "library_prep"
covariates = ["library_prep", "sex", "donor", "region"]
labels_keys = ["cross_area_subclass", "cross_area_cluster"]

cell_count = cell_count.loc[[i in split_value for i in cell_count[split_key]], :]

# 构建组成 AnnData：行 = 样本（library_prep），列 = 细胞类型计数，obs 为协变量
adata = ad.AnnData(sp.sparse.csr_matrix(np.zeros((cell_count.shape[0], 0)), dtype=np.float32))
adata.obs_names = cell_count.index
adata.obs = cell_count

abundances = dat.from_scanpy(
    adata,
    cell_type_identifier=labels_keys[0],
    sample_identifier=random_effect,
)
tmp = adata.obs.loc[:, covariates].drop_duplicates()
tmp.index = tmp.library_prep
tmp = tmp.drop([random_effect], axis=1)
abundances.obs = tmp.loc[abundances.obs.index, :]

# 参考类型逐一轮换：HMC 采样 + 接受率复核，汇总稳健效应
formula = "sex + "
tests = ["region"]
results_table = None
for ct in abundances.var.index:
    model = mod.CompositionalAnalysis(abundances, formula=formula + tests[0], reference_cell_type=ct)
    results = model.sample_hmc()

    accepted = results.sample_stats["is_accepted"].to_numpy()
    if accepted.sum() / accepted.shape[1] < 0.6:
        continue

    current_results = results.summary_prepare()[1].reset_index()
    current_results["Reference Cell Type"] = ct
    results_table = current_results if results_table is None else pd.concat([results_table, current_results], axis=0)

results_table.to_csv("results/sccoda_effects_summary.csv", index=False)
with open("results/sccoda_last_model.pkl", "wb") as f:
    pickle.dump(results, f)
