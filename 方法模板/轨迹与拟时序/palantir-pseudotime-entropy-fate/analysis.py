## Palantir 扩散多尺度流形拟时序、分化熵（Entropy）与终端命运概率推断
import palantir
import scanpy as sc

adata = sc.read_h5ad("Data/scrna_preprocessed.h5ad")
start_cell = adata.obs_names[adata.obs["cell_type"] == "HSC"][0]

# 基于各向异性扩散核计算特征向量，并在特征谱间隙（Eigenspace gap）自适应确定多尺度扩散空间
dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(adata)
imputed = palantir.utils.run_magic_imputation(adata)

# 从指定起始干细胞出发计算各细胞拟时序位置、Shannon 分化熵及向各终末状态转移概率
pr_res = palantir.core.run_palantir(
    adata,
    early_cell=start_cell,
    num_waypoints=1200,
)

adata.obs["palantir_pseudotime"] = pr_res.pseudotime
adata.obs["palantir_entropy"] = pr_res.entropy
adata.obsm["palantir_branch_probs"] = pr_res.branch_probs.to_numpy()
adata.write_h5ad("results/palantir_trajectory.h5ad", compression="gzip")
pr_res.branch_probs.to_csv("results/palantir_branch_probabilities.csv")
