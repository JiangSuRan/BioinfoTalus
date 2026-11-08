## scirpy 单细胞 TCR 链质控、克隆型定义与克隆扩增分析
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir

list_samples = ["STR4", "STR5", "STR21"]
list_basedir = ["Data/STR4/", "Data/STR5/", "Data/STR21/"]

# 逐样本合并表达与 10x VDJ contig 注释
list_adata = []
for d, s in zip(list_basedir, list_samples):
    a = sc.read_10x_mtx(f"{d}{s}_T_5DE/filtered_feature_bc_matrix", cache=True)
    a.obs["sample"] = s
    a.var_names_make_unique()
    f_tcr = f"{d}{s}_T_5DE/filtered_contig_annotations.csv"
    a_tcr = ir.io.read_10x_vdj(f_tcr)
    ir.pp.merge_with_ir(a, a_tcr)
    list_adata.append(a)

adata = list_adata[0].concatenate(list_adata[1:], batch_key="sample")

# 基本预处理
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# 链质控：链配对亚型与受体类型统计
ir.tl.chain_qc(adata)

# 核苷酸精确克隆型（双臂一致、取主链）
ir.pp.ir_dist(adata)
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")

# 氨基酸对齐距离的克隆型簇（宽容定义，跨克隆共享可检出）
ir.pp.ir_dist(adata, metric="alignment", sequence="aa", cutoff=15)
ir.tl.define_clonotype_clusters(adata, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any")

# 克隆扩增与可视化
ir.tl.clonal_expansion(adata)
adata.obs["clone_id_size"] = adata.obs["clone_id_size"].astype("category")
sc.pl.umap(adata, color=["clonal_expansion", "clone_id_size", "chain_pairing"], show=False, save="_clonotype_qc.pdf")

ir.pl.clonal_expansion(adata, groupby="sample", clip_at=4, normalize=False, figsize=(6, 4), show=False)
import matplotlib.pyplot as plt
plt.savefig("results/scirpy_clonal_expansion.pdf", dpi=300, bbox_inches="tight")

adata.write_h5ad("results/scirpy_tcr_annotated.h5ad")
