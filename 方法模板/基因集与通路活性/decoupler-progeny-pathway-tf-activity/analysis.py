## decoupler CollecTRI/PROGENy 网络的 ULM 转录因子与通路活性推断
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc

adata = sc.read_h5ad("Data/fibroblast_annotated.h5ad")

# CollecTRI 转录因子调控网络（无复合物拆分）
net = dc.get_collectri(organism="human", split_complexes=False)

# ULM（univariate linear model）推断每细胞 TF 活性（data 层需已归一化）
dc.run_ulm(
    mat=adata,
    net=net,
    source="source",
    target="target",
    weight="weight",
    verbose=True,
)

acts = dc.get_acts(adata, obsm_key="ulm_estimate")

# 供 PROGENy 通路活性时同样调用 dc.run_mlm（向量模型适配矩阵权重）
progeny_net = dc.get_progeny(organism="human", top=500)
dc.run_mlm(
    mat=adata,
    net=progeny_net,
    source="source",
    target="target",
    weight="weight",
    verbose=True,
)
pathway_acts = dc.get_acts(adata, obsm_key="mlm_estimate")

acts.write_h5ad("results/decoupler_tf_activity.h5ad")
pathway_acts.write_h5ad("results/decoupler_progeny_activity.h5ad")

# 按细胞类型汇总平均 TF 活性
tf_mean = pd.DataFrame(acts.X, index=acts.obs_names, columns=acts.var_names)
tf_mean["celltype"] = acts.obs["celltype"].values
tf_summary = tf_mean.groupby("celltype").mean()
tf_summary.to_csv("results/decoupler_tf_activity_by_celltype.csv")
