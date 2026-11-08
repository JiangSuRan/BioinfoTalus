## pyVIA 惰性随机游走拓扑轨迹、终端状态识别与马尔可夫拟时序
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import pyVIA.core as via

adata = sc.read_h5ad("Data/via_input.h5ad")
df_ids = pd.read_csv("Data/via_ids.csv")
true_label = df_ids["group_id"].tolist()

ncomps, knn, random_seed = 30, 10, 2

# VIA 在低维表征（如 X_css/X_diffusion）上构建 kNN 图并做拓扑图剪枝
v0 = via.VIA(
    adata.obsm["X_css"][:, 0:ncomps],
    true_label,
    jac_std_global=0.15,
    dist_std_local=1,
    knn=knn,
    cluster_graph_pruning_std=1,
    too_big_factor=0.3,
    preserve_disconnected=True,
    dataset="group",
    random_seed=random_seed,
)
v0.run_VIA()

# 流线图 + 马尔可夫拟时序
embedding = adata.obsm["X_umap"]
fig, ax = plt.subplots(figsize=(6, 6))
via.via_streamplot(v0, embedding, scatter_size=5, ax=ax)
plt.axis("equal")
fig.savefig("results/via_streamplot.png", dpi=300, bbox_inches="tight")

pt = v0.single_cell_pt_markov
np.savetxt("results/via_pseudotime.csv", pt, delimiter=",")

adata.obs["via_pseudotime"] = pt
adata.write_h5ad("results/via_trajectory.h5ad")
