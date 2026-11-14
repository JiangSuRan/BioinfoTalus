## dbDiffusion Perturb-seq 效应矩阵构建、扰动/基因共聚类与未见扰动交叉映射
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

adata = sc.read_h5ad("Data/perturb_lognorm.h5ad")
pert_col = "perturbation"
control_label = "non-targeting"

# 从 log-归一化表达重构扰动均值偏移矩阵并截断极端值
X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
labels = adata.obs[pert_col].astype(str).to_numpy()
ctrl_mean = X[labels == control_label].mean(axis=0)

ko_rows, pert_names = [], []
for p in pd.unique(labels):
    if p == control_label:
        continue
    mask = labels == p
    if mask.sum() < 10:
        continue
    ko_rows.append(np.clip(X[mask].mean(axis=0) - ctrl_mean, -1.0, 1.0))
    pert_names.append(p)

KO = np.asarray(ko_rows, dtype=np.float32)

# 扰动级 PCA + Leiden + UMAP（对应原论文 Fig 3a）
pa = sc.AnnData(KO.copy())
pa.obs["pert"] = pd.Series(pert_names, index=pa.obs_names)
sc.pp.scale(pa)
sc.tl.pca(pa, n_comps=min(64, min(KO.shape) - 1), svd_solver="arpack", random_state=0)
pert_emb = pa.obsm["X_pca"]

pe = sc.AnnData(pert_emb.copy())
pe.obs["pert"] = pd.Series(pert_names, index=pe.obs_names)
sc.pp.neighbors(pe, n_neighbors=min(10, len(pert_names) - 1), use_rep="X", random_state=0)
sc.tl.leiden(pe, resolution=0.5, flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(pe, random_state=0)
pert_leiden = pe.obs["leiden"].astype(str)
pert_umap = pe.obsm["X_umap"]

# 基因共表达聚类：载荷内积近似基因相关矩阵后降维聚类（对应 Fig 3c）
ag = adata.copy()
sc.pp.scale(ag)
sc.tl.pca(ag, n_comps=min(64, min(ag.shape) - 1), svd_solver="arpack", random_state=0)
load_mat = ag.varm["PCs"]
k_pc = min(50, load_mat.shape[1])
cor_mat = load_mat[:, :k_pc] @ load_mat[:, :k_pc].T

gc_ad = sc.AnnData(cor_mat.astype(np.float32))
sc.tl.pca(gc_ad, n_comps=min(50, cor_mat.shape[0] - 1), svd_solver="arpack", random_state=0)
gn = sc.AnnData(gc_ad.obsm["X_pca"].copy())
gn.obs_names = adata.var_names.tolist()
sc.pp.neighbors(gn, n_neighbors=min(19, gn.n_obs - 1), random_state=0)
sc.tl.leiden(gn, resolution=0.7, flavor="igraph", n_iterations=2, directed=False, random_state=0)
gene_clusters = gn.obs["leiden"].astype(str)

# 两步交叉映射估计目标/未见扰动嵌入（对应 Fig 8）：定位基因簇 -> 匹配扰动簇质心
target_pert = pert_names[0]
cluster_dict = {
    cl: [pert_names[i] for i, c in enumerate(pert_leiden) if c == cl]
    for cl in sorted(pert_leiden.unique())
}
emb_dict = {pert_names[i]: pert_emb[i] for i in range(len(pert_names))}

if target_pert in gene_clusters.index:
    target_gc = gene_clusters.loc[target_pert]
else:
    pi = pert_names.index(target_pert)
    top_genes = [adata.var_names[i] for i in np.argsort(np.abs(KO[pi]))[::-1] if adata.var_names[i] in gene_clusters.index]
    target_gc = gene_clusters.loc[top_genes[0]]

same_gc = set(gene_clusters[gene_clusters == target_gc].index) - {target_pert}
score_rows = []
for cl, members in cluster_dict.items():
    matched = sum(1 for m in members if m in same_gc)
    prop = matched / len(members)
    score_rows.append({
        "pert_cluster": cl,
        "matched": matched,
        "total": len(members),
        "proportion": prop,
        "score": matched * prop,
    })

score_tbl = pd.DataFrame(score_rows).sort_values("score", ascending=False).reset_index(drop=True)
best_cluster = score_tbl.iloc[0]["pert_cluster"]
est_emb = np.mean([emb_dict[m] for m in cluster_dict[best_cluster]], axis=0)

np.savez(
    "results/dbdiffusion_pert_embeddings.npz",
    pert_names=np.asarray(pert_names),
    pert_emb=pert_emb,
    est_target_emb=est_emb,
    target_pert=np.asarray([target_pert]),
)
pd.DataFrame({
    "perturbation": pert_names,
    "pert_cluster": pert_leiden.values,
    "umap1": pert_umap[:, 0],
    "umap2": pert_umap[:, 1],
}).to_csv("results/dbdiffusion_cluster_assignments.csv", index=False)
score_tbl.to_csv("results/dbdiffusion_crossmap_scores.csv", index=False)
