## Hotspot 细胞状态流形空间自相关基因识别与局部协同调控模块聚类
import hotspot
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("Data/scrna_counts_with_pca.h5ad")

# 基于原始离散 UMI 计数（danb 负二项模型）与 PCA 潜在流形建立 kNN 图并检验局部空间自相关显著基因
hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model="danb",
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts",
)
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations(jobs=1)

hs_genes = hs_results.loc[hs_results["FDR"] < 0.05].index
local_cor = hs.compute_local_correlations(hs_genes, jobs=1)
modules = hs.create_modules(min_gene_threshold=20, core_only=True, fdr_threshold=0.05)
module_scores = hs.calculate_module_scores()

adata.obsm["hotspot_module_scores"] = module_scores.to_numpy()
hs_results.to_csv("results/hotspot_autocorrelation_genes.csv")
pd.DataFrame({"gene": modules.index, "module": modules.values}).to_csv("results/hotspot_gene_modules.csv", index=False)
module_scores.to_csv("results/hotspot_cell_module_scores.csv")
