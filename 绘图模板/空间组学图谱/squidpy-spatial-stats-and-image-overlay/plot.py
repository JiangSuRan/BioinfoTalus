## Squidpy 空间聚落分布、邻域富集、距离共现与 Ripley L 空间聚集联图
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad("Data/spatial_squidpy_stats.h5ad")

sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
sq.gr.nhood_enrichment(adata, cluster_key="celltype", n_jobs=1, show_progress_bar=False)
sq.gr.centrality_scores(adata, cluster_key="celltype", n_jobs=1, show_progress_bar=False)
sq.gr.co_occurrence(adata, cluster_key="celltype", n_jobs=1, show_progress_bar=False)
sq.gr.ripley(adata, cluster_key="celltype", mode="L")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
sq.pl.spatial_scatter(
    adata,
    color="celltype",
    shape=None,
    size=25,
    title="A. Spatial Cell Type Map",
    ax=axes[0],
)
sq.pl.ripley(
    adata,
    cluster_key="celltype",
    mode="L",
    ax=axes[1],
)
axes[1].set_title("B. Ripley's L Spatial Clustering")
fig.subplots_adjust(wspace=0.45)
fig.savefig("results/squidpy_spatial_statistics_suite.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

# nhood_enrichment 内部使用 make_axes_locatable 与 imshow；以 300 dpi PNG 独立落盘可避免多面板裁切及部分 PDF 后端 imshow 失真
sq.pl.nhood_enrichment(
    adata,
    cluster_key="celltype",
    title="Neighborhood Enrichment Z-score",
    figsize=(5.2, 4.6),
)
plt.gcf().savefig("results/squidpy_nhood_enrichment.png", dpi=300, bbox_inches="tight")
plt.close()

# Squidpy 的 co_occurrence 与 centrality_scores 内部自行管理多面板坐标轴，绘图后显式调用 savefig 落盘至 results/
sq.pl.co_occurrence(
    adata,
    cluster_key="celltype",
    clusters=adata.obs["celltype"].cat.categories.tolist(),
    figsize=(12, 4),
)
plt.gcf().savefig("results/squidpy_co_occurrence.pdf", dpi=300, bbox_inches="tight")
plt.close()

sq.pl.centrality_scores(
    adata,
    cluster_key="celltype",
    figsize=(10, 4.5),
)
plt.gcf().savefig("results/squidpy_centrality_scores.pdf", dpi=300, bbox_inches="tight")
plt.close()
