## CellTypist 预训练参考模型单细胞自动化注释与 Majority Voting 聚类校正
import celltypist
import scanpy as sc
from celltypist import models

adata = sc.read_h5ad("Data/scrna_counts.h5ad")

# CellTypist 要求输入经 1e4 深度标准化并 log1p 变换的对数表达矩阵
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

model = models.Model.load(model="Data/Immune_All_Low.pkl")
# 开启 majority_voting 结合局部细胞邻域共表达共识投票，平抑单细胞尺度噪声错分
predictions = celltypist.annotate(
    adata,
    model=model,
    majority_voting=True,
    over_clustering="leiden",
)
adata = predictions.to_adata()

adata.write_h5ad("results/celltypist_annotated.h5ad", compression="gzip")
adata.obs[["predicted_labels", "majority_voting", "conf_score"]].to_csv("results/celltypist_predictions.csv")
