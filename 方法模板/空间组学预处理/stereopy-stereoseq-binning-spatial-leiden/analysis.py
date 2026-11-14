## Stereopy Stereo-seq 芯片 bin 分箱加载、空间质控、空间邻域聚类与标志基因识别
import stereo as st

# 按指定 bin 尺度（如 bin50 约等于单细胞级视野）加载 GEF 空间芯片矩阵并保留原始计数检查点
data = st.io.read_gef(file_path="Data/stereoseq_tissue.gef", bin_size=50)
data.tl.cal_qc()
data.tl.filter_cells(min_gene=200, min_n_genes_by_counts=3, pct_counts_mt=15, inplace=True)
data.tl.raw_checkpoint()

data.tl.normalize_total(target_sum=1e4)
data.tl.log1p()
data.tl.highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5, res_key="highly_variable_genes")
data.tl.scale(max_value=10, zero_center=True)
data.tl.pca(use_highly_genes=True, n_pcs=30, res_key="pca")
data.tl.neighbors(pca_res_key="pca", n_pcs=30, res_key="neighbors")
data.tl.spatial_neighbors(neighbors_res_key="neighbors", res_key="spatial_neighbors")
data.tl.leiden(neighbors_res_key="neighbors", resolution=0.6, res_key="leiden")
data.tl.find_marker_genes(cluster_res_key="leiden", method="wilcoxon_test", res_key="marker_genes")

st.io.write_h5ad(data, use_raw=True, use_result=True, output="results/stereopy_analyzed.h5ad")
