## rapids-singlecell 与 Scanpy 百万级单细胞 GPU 加速质控、黑名单高变基因过滤与 Harmony 整合
import cupy as cp
import numpy as np
import rmm
import scanpy as sc
import rapids_singlecell as rsc
from rmm.allocators.cupy import rmm_cupy_allocator

np.random.seed(66)

# 托管内存允许 GPU 超额分配，适配百万级矩阵
rmm.reinitialize(managed_memory=True)
cp.cuda.set_allocator(rmm_cupy_allocator)

adata = sc.read_h5ad("Data/cima_raw.h5ad")

# 质控协变量：线粒体/核糖体/非编码 RNA/血红蛋白
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
adata.var["rp"] = adata.var_names.str.match(r"^RP[SL][0-9]")
adata.var["ncRNA"] = adata.var_names.str.match(r"^[A-Z][A-Z][0-9].*\.[0-9]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp", "ncRNA", "hb"], percent_top=None, log1p=False, inplace=True)

# 双细胞标记（上游 Scrublet 预测列）与阈值过滤
adata = adata[adata.obs["predicted_doublets"] == False, :]
adata = adata[adata.obs["pct_counts_mt"] < 15, :]
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, max_genes=6000)
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, max_counts=25000)

# HVG 后剔除线粒体/核糖体/非编码/LINC，保留纯净生物学特征
sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="seurat_v3")
remove = adata.var_names.str.match(r"^MT-") | adata.var_names.str.match(r"^RP[SL]") | \
    adata.var_names.str.match(r"^[A-Z][A-Z][0-9].*\.[0-9]") | adata.var_names.str.match(r"(^LOC|LINC)[1-9]*")
adata.var["highly_variable"] = adata.var["highly_variable"] & ~remove

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG 子集缩放 PCA，主成分回写全对象
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()

# GPU 加速 Harmony 按样本批次整合后聚类
rsc.tl.harmony_integrate(adata, "sample")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=19, use_rep="X_pca_harmony")
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=1.5, n_iterations=2, key_added="leiden")

adata.write_h5ad("results/cima_gpu_integrated.h5ad")
