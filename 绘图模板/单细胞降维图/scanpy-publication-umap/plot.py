## Scanpy 出版级高密度单细胞 UMAP 嵌入散点图
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad("data/umap_input.h5ad")

# 高密度点云用小点径与适度透明度，去掉边框与网格保持出版级版面
fig, ax = plt.subplots(figsize=(6, 6))
sc.pl.embedding(
    adata,
    basis="umap",
    color="cell_type",
    size=2,
    alpha=0.6,
    frameon=False,
    title="",
    legend_fontsize=6,
    ax=ax,
    show=False,
)
fig.savefig("result/publication_umap.pdf", dpi=300, bbox_inches="tight")
