## SCENIC+ 与 pycisTarget 单细胞多组学（scRNA + scATAC）增强子驱动调控网络（eRegulon）推断
import pickle

import numpy as np
import scanpy as sc
from pycistarget.motif_enrichment_dem import *
from scenicplus.scenicplus_class import create_SCENICPLUS_object
from scenicplus.wrappers.run_scenicplus import run_scenicplus

adata = sc.read_h5ad("Data/rna_metacell.h5ad")
cistopic_obj = pickle.load(open("Data/cistopic_obj.pkl", "rb"))
menr = pickle.load(open("Data/motifs_menr.pkl", "rb"))

# 剔除非编码 RNA 与 LINC 基因，避免长链非编码基因干扰共表达模块与 TF 匹配
ncRNA_genes = list(adata.var_names[adata.var_names.str.match(r"^[A-Z][A-Z][0-9].*\.[0-9]")])
LINC_genes = list(adata.var_names[adata.var_names.str.match(r"(^LOC|LINC)[1-9]*")])
genes_keep = [g for g in adata.var_names if g not in ncRNA_genes + LINC_genes]
adata = adata[:, genes_keep].copy()

scplus_obj = create_SCENICPLUS_object(
    GEX_anndata=adata,
    cisTopic_obj=cistopic_obj,
    menr=menr,
    ACC_prefix="ACC_",
    GEX_prefix="GEX_",
    bc_transform_func=lambda x: x,
    normalize_imputed_acc=False,
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

# TF 搜索窗为 TSS 上下游 1kb-150kb，eRegulon 需 motif 与共表达双向支持
run_scenicplus(
    scplus_obj=scplus_obj,
    variable=["GEX_cell_type"],
    species="hsapiens",
    assembly="hg38",
    tf_file="Data/utoronto_human_tfs_v_1.01.txt",
    save_path="results/scenicplus_out",
    biomart_host="http://sep2019.archive.ensembl.org/",
    upstream=[1000, 150000],
    downstream=[1000, 150000],
    calculate_TF_eGRN_correlation=False,
    calculate_DEGs_DARs=False,
    export_to_loom_file=False,
    export_to_UCSC_file=False,
    n_cpu=32,
    _temp_dir="Data/temp",
)

pickle.dump(scplus_obj, open("results/scplus_obj.pkl", "wb"))
scplus_obj.uns["eRegulon_metadata"].to_csv("results/eRegulon_metadata.csv")
