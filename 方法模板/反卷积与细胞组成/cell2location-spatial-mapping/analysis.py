## cell2location 单细胞参考回归模型估计与空间贝叶斯反卷积
import numpy as np
import scanpy as sc
import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

# 阶段 1：单细胞参考特征估计（基因过滤 -> 回归模型 -> 后验导出）
adata_sc = sc.read_h5ad("Data/total.singlecell.h5ad")

selected = filter_genes(
    adata_sc,
    cell_count_cutoff=5,
    cell_percentage_cutoff2=0.03,
    nonz_mean_cutoff=1.12,
)
adata_sc = adata_sc[:, selected].copy()

RegressionModel.setup_anndata(adata=adata_sc, batch_key="Samples", labels_key="cell_type")

mod = RegressionModel(adata_sc)
mod.train(max_epochs=250, train_size=1)

adata_sc = mod.export_posterior(
    adata_sc,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500},
)
adata_sc = mod.export_posterior(
    adata_sc,
    use_quantiles=True,
    add_to_varm=["q05", "q50", "q95", "q0001"],
    sample_kwargs={"batch_size": 2500},
)
mod.save("results/reference_signatures/", overwrite=True)
adata_sc.write("results/sc_cell2location.h5ad")

# 阶段 2：空间反卷积（线粒体基因单独建模后剔除 -> 基因对齐 -> 细胞丰度后验）
adata_vis = sc.read_h5ad("Data/st_sample.h5ad")

adata_vis.var["MT_gene"] = [gene.startswith("mt-") for gene in adata_vis.var.index]
adata_vis.obsm["MT"] = adata_vis[:, adata_vis.var["MT_gene"].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]

inf_aver = adata_sc.varm["means_per_cluster_mu_fg"][
    [f"means_per_cluster_mu_fg_{i}" for i in adata_sc.uns["mod"]["factor_names"]]
].copy()
inf_aver.columns = adata_sc.uns["mod"]["factor_names"]

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

cell2location.models.Cell2location.setup_anndata(adata=adata_vis, layer="raw_count")
mod_run = cell2location.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver,
    N_cells_per_location=30,
    detection_alpha=20,
)
mod_run.train(max_epochs=20000, batch_size=None, train_size=1)

adata_vis = mod_run.export_posterior(
    adata_vis,
    sample_kwargs={"num_samples": 1000, "batch_size": mod_run.adata.n_obs},
)

adata_vis.obs[adata_vis.uns["mod"]["factor_names"]] = adata_vis.obsm["q05_cell_abundance_w_sf"]
mod_run.save("results/cell2location_map/", overwrite=True)
adata_vis.write("results/st_cell2location_res.h5ad")
