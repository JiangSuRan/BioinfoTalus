## inferCNVpy 基因组滑动窗口拷贝数变异推断与非整倍体恶变亚群聚类
import infercnvpy as cnv
import scanpy as sc

adata = sc.read_h5ad("Data/scrna_with_genomic_pos.h5ad")

# 利用已知二倍体免疫/基质细胞群作参考基线，沿染色体顺序按 250 基因滑动窗口估算 CNV 状态
cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=["T_cells", "Fibroblasts"],
    window_size=250,
    step=1,
    n_jobs=1,
)
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata, resolution=0.5)
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

adata.write_h5ad("results/infercnvpy_scored.h5ad", compression="gzip")
adata.obs[["cell_type", "cnv_leiden", "cnv_score"]].to_csv("results/infercnvpy_cell_scores.csv")
