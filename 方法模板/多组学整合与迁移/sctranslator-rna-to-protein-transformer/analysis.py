## scTranslator Performer 转录组到蛋白丰度跨模态翻译推理
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn.model_selection import train_test_split

import sys
sys.path.append("code/model")
from performer_enc_dec import PerformerEncDec
from utils import SCDataset, fix_SCDataset, test

# 加载 stage2 预训练检查点
device = "cuda" if torch.cuda.is_available() else "cpu"
pretrain_checkpoint = "Data/checkpoint/stage2_single-cell_scTranslator.pt"
model = torch.load(pretrain_checkpoint, map_location=torch.device(device))

# 测试数据：配对 RNA（2 万基因）与 CITE-seq 蛋白
enc_max_seq_len = 20000
dec_max_seq_len = 1000
test_batch_size = 4

scRNA_adata = sc.read_h5ad("Data/GSM5008737_RNA_finetune_withcelltype.h5ad")[:100]
scP_adata = sc.read_h5ad("Data/GSM5008738_protein_finetune_withcelltype.h5ad")[:100]

test_rna = scRNA_adata
test_protein = scP_adata[test_rna.obs.index]

# 对齐（fix）数据集构建与推理
my_testset = fix_SCDataset(test_rna, test_protein, enc_max_seq_len, dec_max_seq_len)
test_loader = torch.utils.data.DataLoader(my_testset, batch_size=test_batch_size, drop_last=True)

test_loss, test_ccc, y_hat, y = test(model, device, test_loader)

y_pred = pd.DataFrame(y_hat, columns=test_protein.var.index.tolist())
y_truth = pd.DataFrame(y, columns=test_protein.var.index.tolist())

y_pred.to_csv("results/sctranslator_y_pred.csv")
y_truth.to_csv("results/sctranslator_y_truth.csv")
pd.DataFrame({"test_loss": [test_loss], "test_ccc": [test_ccc]}).to_csv("results/sctranslator_performance.csv")
