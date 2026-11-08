## SpatialData 跨平台空间组学统一容器构建、表格提取与导出
import spatialdata
import spatialdata_io
import spatialdata_plot
import scanpy as sc

# xenium reader 将输出目录解析为 Shapes（细胞多边形）/Points（转录本）/Images/Tables 统一结构
sdata = spatialdata_io.xenium("Data/xenium_outs/")
sdata

# 取出单细胞表格并保留计数副本
adata = sdata.tables["table"].copy()
adata.layers["counts"] = adata.X.copy()

# 空间坐标写入 obs 便于与外部工具交换
adata.obs["x"] = adata.obsm["spatial"][:, 0]
adata.obs["y"] = adata.obsm["spatial"][:, 1]

# 计数矩阵（mtx）与元数据导出，供 Seurat 等下游生态使用
from scipy.io import mmwrite
mmwrite("results/xenium_matrix.mtx", adata.layers["counts"])
adata.obs.to_csv("results/xenium_metadata.csv")

with open("results/xenium_gene_names.txt", "w") as f:
    f.writelines([i + "\n" for i in list(adata.var_names)])

sdata.write("results/xenium_spatialdata.zarr")
