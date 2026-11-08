## CellRank 速度核与表达相似性核组合的转移矩阵、命运概率与驱动基因
import cellrank as cr
import scanpy as sc
import scvelo as scv

adata = sc.read_h5ad("Data/velocity_annotated.h5ad")

# 剪接动力学速率与一/二阶矩（CellRank 速度核前置）
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")

# 速度核（动力学方向）与连通性核（表达相似）加权组合为转移矩阵
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck

# 速率场投影（含终态吸收计算由 estimator 完成）
vk.plot_projection(basis="umap", color="celltype", recompute=True, stream=True,
                   save="results/cellrank_velocity_projection.svg")

# 宏观状态与吸收态命运概率
g = cr.estimators.GPCCA(combined_kernel)
g.compute_schur(n_components=10)
g.compute_macrostates(n_states=3, cluster_key="celltype")
g.set_terminal_states_from_macrostates(names=["Terminal_1", "Terminal_2", "Terminal_3"])
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False, save="results/cellrank_fate_probs.png")

fate_probs = g.fate_probabilities
fate_probs.to_adata()
adata.write_h5ad("results/cellrank_fates.h5ad")
