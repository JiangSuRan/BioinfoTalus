## CellRank 终末命运环状单纯形投影、谱系基因趋势 GAM 曲线与级联热图
import cellrank as cr
import matplotlib.pyplot as plt
import scanpy as sc

sc.settings.figdir = "results"
adata = sc.read_h5ad("Data/cellrank_fitted.h5ad")
model = cr.models.GAM(adata)
driver_genes = adata.var_names[:6].tolist()
lineages = adata.obsm["lineages_fwd"].names.tolist()

# 环状投影直观呈现多谱系吸收概率单纯形分布，级联热图沿特定终末谱系潜时间排序驱动基因激活波峰
cr.pl.circular_projection(
    adata,
    keys="clusters",
    legend_loc="right margin",
    save="_cellrank_circular.pdf",
)
cr.pl.gene_trends(
    adata,
    model=model,
    data_key="X",
    genes=driver_genes[:4],
    same_plot=True,
    hide_cells=True,
    save="_cellrank_gene_trends.pdf",
)
cr.pl.heatmap(
    adata,
    model=model,
    genes=driver_genes,
    lineages=lineages[0],
    time_key="latent_time",
    show_all_genes=True,
    save="_cellrank_cascade_heatmap.pdf",
)
