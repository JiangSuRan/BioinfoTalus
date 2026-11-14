## scVelo 四联 RNA 速率可视化（剪切比例/UMAP 速率流线/相图与 PAGA 速度摘要图）
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv

adata = scv.read("Data/scvelo_dynamical_fitted.h5ad")
cluster_key = "clusters"
driver_gene = adata.var["fit_likelihood"].idxmax()

# 展示各亚群未剪切/剪切比例、流线矢量场、动态相图及由速度方向定向的 PAGA 拓扑转移
fig, axes = plt.subplots(2, 2, figsize=(13, 11))
scv.pl.proportions(adata, groupby=cluster_key, ax=axes[0, 0], show=False)
scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color=cluster_key,
    legend_loc="right margin",
    ax=axes[0, 1],
    show=False,
)
scv.pl.scatter(
    adata,
    var_names=driver_gene,
    color=cluster_key,
    ax=axes[1, 0],
    show=False,
)
scv.pl.paga(
    adata,
    basis="umap",
    size=50,
    alpha=0.15,
    min_edge_width=2,
    node_size_scale=1.5,
    ax=axes[1, 1],
    show=False,
)
fig.tight_layout()
fig.savefig("results/scvelo_velocity_paga_suite.pdf", dpi=300, bbox_inches="tight")
