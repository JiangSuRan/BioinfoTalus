## pyVIPER 基于 NaRnEA 与 aREA 调控子富集模型从单细胞表达谱推断主调节蛋白活性
import anndata
import pandas as pd
import pyviper
import scanpy as sc

adata = sc.read_h5ad("Data/scrna_preprocessed.h5ad")
regulon = pyviper.io.read_aracne("Data/aracne_network.tsv")

# 过滤调控子中不在当前单细胞表达空间内的靶基因，避免富集统计量分母偏倚
regulon.filter_targets(adata.var_names)

# 将细胞表达谱转换为基因特征签名并用非参数 NaRnEA 估算蛋白/转录因子调控活性（无背景置换开销）
ges = pyviper.pp.ges_signature(adata, method="zscore")
activity_narnea = pyviper.viper(
    gex_data=ges,
    interactome=regulon,
    enrichment="narnea",
    eset_filter=False,
)

# 按细胞亚群对单细胞蛋白活性执行 Stouffer Z 检验聚合以识别群特异性主调节因子
cluster_meta = pyviper.pp.stouffer(
    activity_narnea,
    cluster_col="cell_type",
)

activity_narnea.write_h5ad("results/pyviper_protein_activity.h5ad", compression="gzip")
cluster_meta.to_df().T.to_csv("results/pyviper_cluster_master_regulators.csv")
