## GenKI VGAE 基因调控网络虚拟敲除与潜空间 KL 距离排序
import numpy as np
import scanpy as sc

import GenKI as gk
from GenKI.dataLoader import DataLoader
from GenKI.preprocessing import build_adata
from GenKI.train import VGAE_trainer
from GenKI import utils

adata = build_adata("Data/microglial_WT.h5ad")
adata = adata[:, :].copy()

# 构建主成分回归 GRN 并准备 WT / 虚拟 KO（TREM2 置零）图数据
data_wrapper = DataLoader(
    adata,
    target_gene=["TREM2"],
    target_cell=None,
    obs_label="ident",
    GRN_file_dir="results/GRNs",
    rebuild_GRN=True,
    n_cpus=8,
)
data_wt = data_wrapper.load_data()
data_ko = data_wrapper.load_kodata()

# VGAE 训练（WT 图）
hyperparams = {"epochs": 100, "lr": 7e-4, "beta": 1e-4, "seed": 8096}
sensei = VGAE_trainer(
    data_wt,
    epochs=hyperparams["epochs"],
    lr=hyperparams["lr"],
    log_dir=None,
    beta=hyperparams["beta"],
    seed=hyperparams["seed"],
    verbose=False,
)
sensei.train()

# WT/KO 潜变量 KL 距离 → 全基因组差异响应基因排序
z_mu_wt, z_std_wt = sensei.get_latent_vars(data_wt)
z_mu_ko, z_std_ko = sensei.get_latent_vars(data_ko)
dis = gk.utils.get_distance(z_mu_ko, z_std_ko, z_mu_wt, z_std_wt, by="KL")

res_raw = utils.get_generank(data_wt, dis, rank=True)
res_raw.to_csv("results/genki_virtual_ko_rank.csv")

# 置换检验 p 值（打乱 KO 图重训的对照分布）
res_perm = utils.permutation_test(
    data_wrapper, sensei, by="KL", nPerm=1000,
    type2="file", save_path="results/genki_permutation.csv",
)
res_perm.to_csv("results/genki_permutation_summary.csv")
