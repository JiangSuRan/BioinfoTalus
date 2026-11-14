## inferCNVpy 染色体拷贝数变异分级热图与 CNV 得分/亚群多面板 UMAP 图
import infercnvpy as cnv
import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad("Data/infercnvpy_scored.h5ad")

# 按 CNV 聚簇汇总展示常染色体区段增益/缺失图谱，区分二倍体参考群与非整倍体肿瘤群
cnv.pl.chromosome_heatmap(
    adata,
    groupby="cnv_leiden",
    show=False,
)
plt.gcf().savefig("results/infercnvpy_chromosome_heatmap.pdf", dpi=300, bbox_inches="tight")
plt.close()

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
cnv.pl.umap(adata, color="cnv_leiden", ax=axes[0], show=False)
cnv.pl.umap(adata, color="cnv_score", cmap="viridis", ax=axes[1], show=False)
sc.pl.umap(adata, color="cell_type", ax=axes[2], show=False)
fig.tight_layout()
fig.savefig("results/infercnvpy_cnv_umap_panel.pdf", dpi=300, bbox_inches="tight")
