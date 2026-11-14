## Scanorama 多切片空间转录组全景批次矫正与 scRNA-seq 余弦距离标签迁移
import numpy as np
import scanpy as sc
import scanorama
from scipy.spatial.distance import cdist

adata_sc = sc.read_h5ad("Data/scrna_reference.h5ad")
adata_st1 = sc.read_h5ad("Data/visium_slice1.h5ad")
adata_st2 = sc.read_h5ad("Data/visium_slice2.h5ad")

# 利用 Panoramic stitching 互近邻配对同时整合多张空间切片与单细胞参考至统一低维语义空间
adatas = [adata_sc, adata_st1, adata_st2]
adatas_cor = scanorama.correct_scanpy(adatas, return_dimred=True)

adata_merged = sc.concat(
    adatas_cor,
    label="batch",
    keys=["scRNA", "Slice1", "Slice2"],
    uns_merge="unique",
)
sc.pp.neighbors(adata_merged, use_rep="X_scanorama")
sc.tl.umap(adata_merged)
sc.tl.leiden(adata_merged, resolution=0.5, key_added="scanorama_leiden")

# 在 Scanorama 对齐嵌入空间中计算空间 Spot 与单细胞参考群的余弦相似度并分配最近邻类型标签
emb_sc = adatas_cor[0].obsm["X_scanorama"]
emb_sp = np.vstack([adatas_cor[1].obsm["X_scanorama"], adatas_cor[2].obsm["X_scanorama"]])
dist_mat = cdist(emb_sp, emb_sc, metric="cosine")
nearest_idx = np.argmin(dist_mat, axis=1)
inferred_labels = adatas_cor[0].obs["cell_type"].iloc[nearest_idx].to_numpy()

sp_merged = sc.concat([adatas_cor[1], adatas_cor[2]], label="slice", keys=["Slice1", "Slice2"])
sp_merged.obs["predicted_celltype"] = inferred_labels

adata_merged.write_h5ad("results/scanorama_integrated_sc_spatial.h5ad", compression="gzip")
sp_merged.write_h5ad("results/scanorama_spatial_transferred.h5ad", compression="gzip")
