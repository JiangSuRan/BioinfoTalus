## 预测驱动推断 (PPI) 生成扰动去偏差均值、诚实标准误与差异基因检验
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm, pearsonr

adata = ad.read_h5ad("Data/perturb_lognorm.h5ad")
pert_key = "perturbation"
ctrl_label = "non-targeting"
target_pert = "TARGET_KO"

X_all = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
labels = adata.obs[pert_key].astype(str).to_numpy()

ctrl_X = X_all[labels == ctrl_label]
ctrl_mean = ctrl_X.mean(axis=0)
ctrl_se = ctrl_X.std(axis=0, ddof=1) / np.sqrt(ctrl_X.shape[0])

obs_X = X_all[labels == target_pert]
obs_mean = obs_X.mean(axis=0)
obs_se = obs_X.std(axis=0, ddof=1) / np.sqrt(obs_X.shape[0])

target_gen = np.load("Data/generated/TARGET_KO_generated.npy")
naive_mean = target_gen.mean(axis=0)
naive_se = target_gen.std(axis=0, ddof=1) / np.sqrt(target_gen.shape[0])

neighbors = pd.read_csv("Data/cluster_neighbors.csv")["perturbation"].astype(str).tolist()
n_genes = adata.n_vars
correction = np.zeros(n_genes)
bias_dict = {}
means_nb, sizes_nb = [], []
pooled = np.zeros(n_genes)

for nb in neighbors:
    nb_obs = X_all[labels == nb]
    nb_gen = np.load(f"Data/generated/{nb}_generated.npy")
    nb_obs_mean = nb_obs.mean(axis=0)
    bias = nb_obs_mean - nb_gen.mean(axis=0)
    bias_dict[nb] = bias
    correction += bias
    means_nb.append(nb_obs_mean)
    n_k = nb_obs.shape[0]
    pooled += nb_obs.var(axis=0, ddof=0) * n_k
    sizes_nb.append(n_k)

K = len(neighbors)
correction /= K
ppi_mean = naive_mean + correction

# 组内池化方差 + 组间均值方差联合估计 PPI 标准误，避免生成样本过窄导致过度自信
s2 = pooled / max(sum(sizes_nb) - K, 1)
t2 = np.var(np.vstack(means_nb), axis=0, ddof=1) if K > 1 else np.zeros(n_genes)
ppi_var = sum(s2 / n_k + t2 for n_k in sizes_nb)
ppi_se = np.sqrt(ppi_var) / K

obs_effect = obs_mean - ctrl_mean
naive_effect = naive_mean - ctrl_mean
ppi_effect = ppi_mean - ctrl_mean
abs_eff = np.abs(obs_effect)
top_idx = np.argsort(abs_eff)[-min(100, n_genes):]

r_naive = pearsonr(naive_effect[top_idx], obs_effect[top_idx])[0]
r_ppi = pearsonr(ppi_effect[top_idx], obs_effect[top_idx])[0]

z_stat = (ppi_mean - ctrl_mean) / np.sqrt(ppi_se**2 + ctrl_se**2 + 1e-12)
p_values = 2.0 * norm.sf(np.abs(z_stat))
bonf_sig = p_values < (0.05 / n_genes)

z_crit = norm.ppf(0.975)
ov_ppi = np.mean(
    ((ppi_mean - z_crit * ppi_se) <= (obs_mean + z_crit * obs_se))
    & ((ppi_mean + z_crit * ppi_se) >= (obs_mean - z_crit * obs_se))
)
ov_naive = np.mean(
    ((naive_mean - z_crit * naive_se) <= (obs_mean + z_crit * obs_se))
    & ((naive_mean + z_crit * naive_se) >= (obs_mean - z_crit * obs_se))
)

pd.DataFrame({
    "gene": adata.var_names.astype(str),
    "ctrl_mean": ctrl_mean,
    "obs_mean": obs_mean,
    "obs_se": obs_se,
    "naive_mean": naive_mean,
    "ppi_mean": ppi_mean,
    "ppi_se": ppi_se,
    "z_stat": z_stat,
    "pval": p_values,
    "sig_bonf": bonf_sig,
    "abs_effect": abs_eff,
}).to_csv("results/ppi_gene_summary.csv", index=False)

pd.DataFrame(bias_dict, index=adata.var_names.astype(str)).to_csv("results/ppi_neighbor_bias.csv")

pd.DataFrame([{
    "target": target_pert,
    "n_neighbors": K,
    "pearson_naive_top": r_naive,
    "pearson_ppi_top": r_ppi,
    "ci95_overlap_naive": ov_naive,
    "ci95_overlap_ppi": ov_ppi,
    "n_sig_bonferroni": int(bonf_sig.sum()),
}]).to_csv("results/ppi_metrics_summary.csv", index=False)
