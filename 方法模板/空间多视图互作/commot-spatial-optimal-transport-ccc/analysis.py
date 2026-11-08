## COMMOT 最优传输距离约束的空间细胞通讯推断与配体受体流向分析
import scanpy as sc
import squidpy as sq
import commot as ct

adata = sq.read.visium("Data/visium_sample/outs")
adata.var_names_make_unique()
adata.raw = adata
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

adata_dis500 = adata.copy()

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.4)

# CellChat 配受体库经表达率过滤后进入最优传输推断（500 单位距离阈值）
df_cellchat = ct.pp.ligand_receptor_database(species="human", signaling_type=None, database="CellChat")
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)

ct.tl.spatial_communication(
    adata_dis500,
    database_name="cellchat",
    df_ligrec=df_cellchat_filtered,
    dis_thr=500,
    heteromeric=True,
    pathway_sum=True,
)

# 每通路的发送方向场与聚类水平网络（置换检验）
pathways_to_plot = ["FGF", "WNT", "NOTCH"]
adata_dis500.obs["leiden"] = adata.obs["leiden"]
for pathway_use in pathways_to_plot:
    ct.tl.communication_direction(adata_dis500, database_name="cellchat", pathway_name=pathway_use, k=5)
    ct.tl.cluster_communication(
        adata_dis500, database_name="cellchat", pathway_name=pathway_use,
        clustering="leiden", n_permutations=100,
    )
ct.tl.cluster_position(adata_dis500, clustering="leiden")

adata_dis500.write_h5ad("results/commot_spatial_ccc.h5ad")
