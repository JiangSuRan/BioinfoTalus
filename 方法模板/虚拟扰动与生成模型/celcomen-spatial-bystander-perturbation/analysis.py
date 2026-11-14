## Celcomen CCE/SCE 特定供体细胞类型虚拟敲除与邻域旁观者细胞空间差异响应
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial.distance import cdist, pdist, squareform
import torch
import torch_geometric
from torch_geometric.loader import DataLoader as PyGDataLoader

from celcomen.models.celcomen import celcomen
from celcomen.models.simcomen import simcomen
from celcomen.training_plan.train import train
from celcomen.utils.helpers import calc_sphex, normalize_g2g

device = "cuda"
adata = sc.read_h5ad("Data/spatial_tumor.h5ad")
X = np.asarray(adata.X.toarray() if sp.issparse(adata.X) else adata.X, dtype=np.float32)
coords = np.asarray(adata.obsm["spatial"])[:, :2].astype(np.float32)

# 每细胞表达按 L2 归一化，使能量模型在基因交互方向而非测序深度尺度上学习胞间/胞内矩阵
x_tensor = torch.from_numpy(X).to(device)
norm_factor = torch.clamp(torch.sqrt(torch.pow(x_tensor, 2).sum(1)).reshape(-1, 1), min=1e-8)
x_norm = torch.div(x_tensor, norm_factor)

dist_mat = squareform(pdist(coords))
edge_index = torch.from_numpy(np.array(np.where((dist_mat < 15.0) & (dist_mat != 0)))).to(device)

n_genes = adata.n_vars
data = torch_geometric.data.Data(
    x=x_norm,
    pos=torch.from_numpy(coords).to(device),
    y=torch.tensor([0], device=device),
    edge_index=edge_index,
)
loader = PyGDataLoader([data], batch_size=1, shuffle=False)

np.random.seed(0)
init_g2g = np.random.uniform(size=(n_genes, n_genes)).astype("float32")
init_g2g = normalize_g2g((init_g2g + init_g2g.T) / 2)

cce = celcomen(input_dim=n_genes, output_dim=n_genes, n_neighbors=6, seed=0)
cce.set_g2g(torch.from_numpy(init_g2g).to(device))
cce.set_g2g_intra(torch.from_numpy(init_g2g.copy()).to(device))
cce.to(device)
train(200, 2e-2, cce, loader, zmft_scalar=1e-1, seed=1, device=device)

# 定位供体细胞类型 A 阳性细胞执行体外基因敲除置零，并按 30um 半径划定旁观者细胞 B 的近端与远端基线
gene_names = adata.var_names.tolist()
perturb_gene = "IFITM3"
gi = gene_names.index(perturb_gene)
ct = adata.obs["cell_type"].astype(str).to_numpy()
mask_A_ko = (ct == "T cell") & (X[:, gi] > 0)
prompt_x = X.copy()
prompt_x[mask_A_ko, gi] = 0.0

dmin = cdist(coords[mask_A_ko], coords).min(axis=0)
is_B = ct == "Tumor/Astro"
B_near = is_B & (dmin > 0.0) & (dmin < 30.0)
B_far = is_B & (dmin >= 30.0)

# 将学习得到的 g2g 矩阵迁移到 SCE 进行能量最小化迭代采样反事实空间转录组
xp = torch.from_numpy(prompt_x).to(device)
xp = torch.div(xp, torch.clamp(torch.sqrt(torch.pow(xp, 2).sum(1)).reshape(-1, 1), min=1e-8))
sce = simcomen(input_dim=n_genes, output_dim=n_genes, n_neighbors=6, seed=0)
sce.set_g2g(cce.conv1.lin.weight.clone().detach())
sce.set_g2g_intra(cce.lin.weight.clone().detach())
sce.set_sphex(calc_sphex(xp))
sce.to(device)

optimizer = torch.optim.SGD(sce.parameters(), lr=1e-3, momentum=0)
sce.train()
torch.manual_seed(0)
sce_init = sce.gex.clone().detach().cpu().numpy()

for _ in range(50):
    msg, msg_intra, log_z_mft = sce(edge_index, 1)
    loss = -(
        -log_z_mft
        + 1e-1 * torch.trace(torch.mm(msg, torch.t(sce.gex)))
        + 1e-1 * torch.trace(torch.mm(msg_intra, torch.t(sce.gex)))
    )
    loss.backward()
    optimizer.step()
    optimizer.zero_grad()

sce_output = sce.gex.clone().detach().cpu().numpy()
diff = sce_output - sce_init

# 提取受体细胞类型 B 的近 A vs 远 A 差异响应，并与不分细胞类型的通用邻域差异对照
de_B = pd.Series(diff[B_near].mean(0) - diff[B_far].mean(0), index=gene_names).sort_values()
near_all = (dmin > 0.0) & (dmin < 30.0) & (~mask_A_ko)
far_all = (dmin >= 30.0) & (~(ct == "T cell"))
de_generic = pd.Series(diff[near_all].mean(0) - diff[far_all].mean(0), index=gene_names)

z_B = (de_B - de_B.mean()) / (de_B.std() + 1e-12)
pd.DataFrame({
    "gene": de_B.index,
    "de_near_minus_far": de_B.values,
    "de_generic": de_generic.reindex(de_B.index).values,
    "z": z_B.values,
}).to_csv("results/celcomen_bystander_de.csv", index=False)

zone = np.array(["other"] * adata.n_obs, dtype=object)
zone[B_far] = "B far"
zone[B_near] = "B near A (<30um)"
zone[mask_A_ko] = "A perturbed"
pd.DataFrame({
    "cell_id": adata.obs_names.astype(str),
    "x": coords[:, 0],
    "y": coords[:, 1],
    "cell_type": ct,
    "is_B": is_B,
    "dist_to_A": dmin,
    "zone": zone,
    "mean_abs_diff": np.abs(diff).mean(axis=1),
}).to_csv("results/celcomen_bystander_zones.csv", index=False)

np.save("results/celcomen_g2g_inter.npy", cce.conv1.lin.weight.detach().cpu().numpy())
np.save("results/celcomen_sce_output.npy", sce_output)
