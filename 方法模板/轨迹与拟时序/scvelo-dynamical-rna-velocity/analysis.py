## scVelo 单细胞剪接动力学 RNA 速率（RNA Velocity）与潜时间推断
import scanpy as sc
import scvelo as scv

adata = sc.read_h5ad("Data/rna_velocity_input.h5ad")

# 动力学模型依赖 spliced/unspliced 计数层；先归一化再估计一/二阶矩
scv.pp.filter_and_normalize(adata, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# 动力学模型联合估计剪接速率 alpha/beta/gamma，再投影为细胞速率向量
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=8)
scv.tl.latent_time(adata)

adata.write_h5ad("results/scvelo_dynamical.h5ad")
