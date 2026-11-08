## SEACells 基于低维流形核矩阵的单细胞 Metacell 聚合与伪多组学构建
import numpy as np
import scanpy as sc
import scipy.sparse
import sklearn.utils.extmath
from SEACells import core
from sklearn.preprocessing import normalize

ad_obj = sc.read_h5ad("Data/sc_input.h5ad")

# TF-IDF + log1p 后随机 SVD 得 50 维 LSI，逐细胞标准化消除测序深度差
X = ad_obj.X
if scipy.sparse.issparse(X):
    idf = X.shape[0] / X.sum(axis=0)
    X_norm = X.multiply(1 / X.sum(axis=1)).multiply(idf)
else:
    idf = X.shape[0] / X.sum(axis=0)
    X_norm = X / X.sum(axis=1, keepdims=True) * idf
X_norm = normalize(X_norm, norm="l1")
X_lsi = sklearn.utils.extmath.randomized_svd(np.log1p(X_norm * 1e4), 50, random_state=666)[0]
X_lsi -= X_lsi.mean(axis=1, keepdims=True)
X_lsi /= X_lsi.std(axis=1, ddof=1, keepdims=True)
ad_obj.obsm["X_lsi"] = X_lsi

# Metacell 数约为细胞数 1/50，核矩阵建立在 LSI 流形上
n_metacells = int(np.floor(ad_obj.n_obs / 50)) + 1
model = core.SEACells(
    ad_obj,
    build_kernel_on="X_lsi",
    n_SEACells=n_metacells,
    n_waypoint_eigs=10,
    convergence_epsilon=1e-5,
    use_gpu=False,
)
model.construct_kernel_matrix()
model.initialize_archetypes()
model.fit(min_iter=10, max_iter=100)

metacell_adata = core.summarize_by_SEACell(ad_obj, SEACells_label="SEACell", summarize_layer="raw")
metacell_adata.write_h5ad("results/metacell_summarized.h5ad")
