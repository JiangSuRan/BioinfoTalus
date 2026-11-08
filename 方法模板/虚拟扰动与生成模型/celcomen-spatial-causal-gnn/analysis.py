## Celcomen 因果图神经网络细胞间/细胞内基因-基因调控学习与空间反事实扰动
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from scipy.spatial.distance import cdist
from sklearn.neighbors import kneighbors_graph

from celcomen.datareaders.datareader import get_dataset_loaders
from celcomen.models.celcomen import celcomen
from celcomen.models.simcomen import simcomen
from celcomen.training_plan.train import train
from celcomen.utils.helpers import calc_sphex, normalize_g2g

# Xenium 数据预处理（HVG 全保留）
adata = sc.read_10x_h5("Data/xenium/cell_feature_matrix.h5")
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_counts=100)
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

df = pd.read_csv("Data/xenium/cells.csv.gz", index_col=0)
df = df.loc[adata.obs.index]
adata.obsm["spatial"] = df.loc[adata.obs.index, ["x_centroid", "y_centroid"]].values
mask = (df["x_centroid"] < 7000) & (df["x_centroid"] > 6700) & (df["y_centroid"] < 8300) & (df["y_centroid"] > 8000)
adata = adata[mask].copy()

adata.uns["log1p"]["base"] = np.e
sc.pp.highly_variable_genes(adata, flavor="seurat", min_mean=0.5, max_mean=7.5, min_disp=0.5, n_top_genes=adata.shape[1])
hvgs = adata.var.index[adata.var["highly_variable"]]
adata = adata[:, hvgs].copy()
adata.obs["sample"] = "0"
adata.write("results/expr.h5ad")

# 空间 kNN dataloader（k=6 近邻，15 单位距离阈值）
sample = "sample"
k = 6
distance = 15
device = "cuda"
dataloader = get_dataset_loaders(
    "results/expr.h5ad", sample_id_name=sample,
    n_neighbors=k, distance=distance, device=device, verbose=False,
)

# Celcomen GNN：细胞间 g2g 与细胞内 g2g_intra 联合学习
n = len(hvgs)
lr = 1e-1
zmft = 1e-1
seed = 0
epochs = 200
model = celcomen(input_dim=n, output_dim=n, n_neighbors=k, seed=seed)
input_g2g = np.random.uniform(size=(n, n)).astype("float32")
input_g2g = normalize_g2g((input_g2g + input_g2g.T) / 2)
model.set_g2g(torch.from_numpy(input_g2g))
model.set_g2g_intra(torch.from_numpy(input_g2g))
model.to(device)

losses = train(epochs, lr, model, dataloader, zmft_scalar=zmft, seed=seed, device=device)

# 空间反事实：单 spot 目标基因置零后以 simcomen 采样邻域响应
expr = torch.from_numpy(adata[:, hvgs].X.toarray() if hasattr(adata[:, hvgs].X, "toarray") else adata[:, hvgs].X)
norm_factor = torch.sqrt(torch.pow(expr, 2).sum(1)).reshape(-1, 1)
expr = torch.div(expr, norm_factor).to(device)

pos = torch.from_numpy(adata.obsm["spatial"])
edges = kneighbors_graph(pos, k, include_self=False).todense()
edges = torch.from_numpy(np.array(np.where(edges))).to(device)

model_sim = simcomen(input_dim=n, output_dim=n, n_neighbors=k, seed=seed)
model_sim.set_g2g(model.conv1.lin.weight.clone().detach())
model_sim.set_g2g_intra(model.lin.weight.clone().detach())
input_sphex = calc_sphex(expr).clone().detach().numpy()
model_sim.set_sphex(torch.from_numpy(input_sphex.astype("float32")))
model_sim.to(device)

fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(losses, lw=2, color="#fe86a4")
ax.set_xlim(0, epochs)
ax.set(xlabel="Epochs", ylabel="Model loss")
fig.savefig("results/celcomen_loss.png", dpi=300, bbox_inches="tight")

torch.save(model.state_dict(), "results/celcomen_model.pt")
