## SIRV 将 scRNA 剪接/未剪接矩阵迁移至空间坐标系估算空间 RNA 速率矢量场
import scanpy as sc
import scvelo as scv
import SIRV

adata_st = sc.read_h5ad("Data/spatial_preprocessed.h5ad")
adata_sc = sc.read_h5ad("Data/scrna_unspliced_spliced.h5ad")

common_genes = list(set(adata_st.var_names).intersection(set(adata_sc.var_names)))
adata_st_sub = adata_st[:, common_genes].copy()
adata_sc_sub = adata_sc[:, common_genes].copy()

# SIRV 基于 SpaGE 领域适应空间投影，将参考单细胞中的 spliced/unspliced 表达插补至空间 Spot
adata_sirv = SIRV.ed.SIRV(adata_st_sub, adata_sc_sub, n_pv=30)
adata_sirv.obsm["X_xy_loc"] = adata_st_sub.obsm["spatial"]

scv.pp.moments(adata_sirv, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_sirv, mode="stochastic")
scv.tl.velocity_graph(adata_sirv)
scv.tl.velocity_embedding(adata_sirv, basis="xy_loc")

adata_sirv.write_h5ad("results/sirv_spatial_velocity.h5ad", compression="gzip")
