## dynamo 随机动力学 RNA 速率、连续向量场重建与细胞水平置信度过滤
import dynamo as dyn
import scvelo as scv
import scanpy as sc

adata = sc.read_h5ad("Data/velocity_annotated.h5ad")

# dynamo 预处理：含细胞周期打分的标准化与特征选择
preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
preprocessor.preprocess_adata(adata)

# 随机动力学模型估计 alpha/beta/gamma 剪接动力学参数
dyn.tl.dynamics(adata, model="stochastic", cores=3)
dyn.tl.reduceDimension(adata)

# Pearson 相关法在降维上估计细胞速率（sqrt 变换抑制长尾）
dyn.tl.cell_velocities(adata, method="pearson", other_kernels_dict={"transform": "sqrt"})
dyn.tl.cell_wise_confidence(adata)

# 按细胞类型分组筛选高置信速率
dyn.tl.confident_cell_velocities(adata, group="celltype")

scv.pl.velocity_embedding_stream(
    adata, basis="umap", color=["celltype"], save="results/dynamo_velocity_stream.svg", title="",
)

adata.write_h5ad("results/dynamo_vectorfield.h5ad")
