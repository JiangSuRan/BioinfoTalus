## Squidpy 空间邻接图构建与邻域富集分析
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad("Data/spatial_input.h5ad")

# generic 坐标类型适配非网格排布的空间平台（Visium / Slide-seq / Stereo-seq 通用）
sq.gr.spatial_neighbors(adata, coord_type="generic")

# 邻域富集 z-score 度量各细胞类型对的空间共定位强度，NaN 置 0 便于下游矩阵化使用
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
zscore = adata.uns["cluster_nhood_enrichment"]["zscore"]
adata.uns["cluster_nhood_enrichment"]["zscore"] = np.nan_to_num(zscore)

# 共定位概率按距离区间估计细胞类型对的邻近发生概率
sq.gr.co_occurrence(adata, cluster_key="cluster")

adata.write_h5ad("results/spatial_nhood.h5ad")
