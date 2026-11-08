## 10x Xenium 原位空间转录组数据读取、质控过滤与空间坐标规范化
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

adata = sc.read_10x_h5("Data/cell_feature_matrix.h5")
adata.obs = pd.read_csv("Data/cells.csv.gz", index_col=0)

# 保留 gene_ids 索引并另存基因符号列，兼容两类基因标识的下游匹配
adata.var = adata.var.reset_index().rename(columns={"index": "gene_name"}).set_index("gene_ids")
adata.var.index.name = None

# Xenium 官方 QC flag 标记的控制探针不在读取时剔除，特征保留状态显式记录
adata.var["feature_is_filtered"] = False

# 细胞质心坐标写入 obsm['spatial']，作为下游空间图构建的标准接口
adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
adata.X = csr_matrix(adata.X)

# 过滤低检出细胞与低检出基因，阈值保守，迁移时按数据分布重评
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=10)

adata.write_h5ad("results/xenium_processed.h5ad")
