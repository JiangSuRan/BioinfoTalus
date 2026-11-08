## DestVI 条件 VAE 单细胞亚群状态建模与空间连续反卷积
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scvi.model import CondSCVI, DestVI
import destvi_utils

# 单细胞参考：counts 层 + 2000 HVG（seurat_v3）+ 归一化
sc_adata = sc.read_h5ad("Data/scRNA_sample.h5ad")
sc.pp.filter_genes(sc_adata, min_counts=10)
G = 2000
sc_adata.layers["counts"] = sc_adata.X.copy()
sc.pp.highly_variable_genes(sc_adata, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3")
sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

# 空间数据：counts 层保留原始计数，坐标转 float
st_adata = sc.read_h5ad("Data/visium_sample.h5ad")
st_adata.layers["counts"] = st_adata.X.copy()
sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata
st_adata.obsm["spatial"] = st_adata.obsm["spatial"].astype("float")

# 两侧基因取交集
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()

# CondSCVI 条件模型（按细胞类型标签）训练后迁移到空间 DestVI
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="Annotation")
sc_model = CondSCVI(sc_adata, weight_obs=False)
sc_model.train()

DestVI.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.train(max_epochs=2500)

# 细胞类型比例与亚群状态 gamma 潜变量导出
st_adata.obsm["proportions"] = st_model.get_proportions()
st_adata.obsm["proportions"].to_csv("results/CellProp_DestVI.csv")

ct_thresholds = destvi_utils.automatic_proportion_threshold(st_adata, kind_threshold="primary")
st_adata.uns["ct_thresholds"] = ct_thresholds

for ct, g in st_model.get_gamma().items():
    st_adata.obsm[f"{ct}_gamma"] = g

st_adata.write_h5ad("results/destvi_spatial.h5ad")
