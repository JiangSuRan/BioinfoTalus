## PHATE 多尺度信息几何降维与 MAGIC 马尔可夫亲和图稀疏基因表达去噪插补
import magic
import phate
import scanpy as sc

adata = sc.read_h5ad("Data/scrna_preprocessed.h5ad")

# MAGIC 通过扩散算子 t 次幂在邻近细胞流形上平滑 Dropout 零膨胀，恢复基因共表达连续梯度
magic_op = magic.MAGIC(knn=5, t="auto", n_jobs=1, random_state=42)
X_magic = magic_op.fit_transform(adata.to_df())
adata.layers["MAGIC_imputed"] = X_magic.to_numpy()

# PHATE 以潜在扩散概率距离构造信息熵散度并 MDS 嵌入，同时保留全局分支轨迹拓扑与局部细节
phate_op = phate.PHATE(n_components=2, knn=5, decay=40, t="auto", n_jobs=1, random_state=42)
adata.obsm["X_phate"] = phate_op.fit_transform(adata.X)

adata.write_h5ad("results/phate_magic_denoised.h5ad", compression="gzip")
