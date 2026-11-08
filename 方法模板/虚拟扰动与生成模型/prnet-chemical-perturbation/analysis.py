## PRnet 化学分子结构条件化的单细胞药物扰动转录响应预测
import numpy as np
import pandas as pd
import scanpy as sc

from trainer.PRnetTrainer import PRnetTrainer

adata = sc.read_h5ad("Data/Sci_Plex.h5ad")
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# PRnetTrainer 封装模型构建与数据划分；载入训练好的检查点对留出化合物/剂量组合预测
trainer = PRnetTrainer(
    adata,
    batch_size=512,
    comb_num=1,
    split_key="random_split_0",
    model_save_dir="results/prnet/",
    results_save_dir="results/prnet/",
    x_dimension=5000,
    hidden_layer_sizes=[128],
    z_dimension=64,
    adaptor_layer_sizes=[128],
    comb_dimension=64,
    drug_dimension=1024,
    dr_rate=0.05,
    n_genes=50,
    loss=["GUSS"],
    obs_key="cov_drug",
)
trainer.test("Data/sciplex_best_epoch_all.pt")

# 预测输出组装为 AnnData：行为留出条件的细胞，列为 landmark 基因
pre_array = np.genfromtxt("results/prnet/random_split_0_y_pre_array.csv", delimiter=",")
cov_drug = [x.strip() for x in open("results/prnet/random_split_0_cov_drug_array.csv")]
pred = sc.AnnData(X=pre_array, obs=pd.DataFrame({"cov_drug": pd.Categorical(cov_drug)}))
pred.write_h5ad("results/prnet_predictions.h5ad")
