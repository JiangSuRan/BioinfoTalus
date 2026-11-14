## scGPT 单细胞基础模型基因 Token 化与批次整合
import json
from pathlib import Path

import numpy as np
import scanpy as sc
import torch
from scipy.sparse import issparse

import scgpt
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

set_seed(42)

model_dir = Path("Data/scGPT_bc")
adata = sc.read_h5ad("Data/pbmc_batched.h5ad")
adata.var["gene_name"] = adata.var.index.tolist()
adata.obs["str_batch"] = adata.obs["batch"].astype(str)
adata.obs["batch_id"] = adata.obs["str_batch"].astype("category").cat.codes.values

# 词表对齐：仅保留预训练词表覆盖的基因
vocab = GeneVocab.from_file(model_dir / "vocab.json")
adata.var["id_in_vocab"] = [1 if g in vocab else -1 for g in adata.var["gene_name"]]
adata = adata[:, adata.var["id_in_vocab"] >= 0]

# 归一化 + log1p + HVG 子集 + 基因表达分箱（51 个 token）
preprocessor = Preprocessor(
    use_key="X",
    filter_gene_by_counts=3,
    normalize_total=1e4,
    result_normed_key="X_normed",
    log1p=True,
    result_log1p_key="X_log1p",
    subset_hvg=1200,
    hvg_flavor="cell_ranger",
    binning=51,
    result_binned_key="X_binned",
)
preprocessor(adata, batch_key="str_batch")

# Token 化：每细胞生成基因 ID 序列与分箱表达值
all_counts = (
    adata.layers["X_binned"].toarray()
    if issparse(adata.layers["X_binned"])
    else adata.layers["X_binned"]
)
genes = adata.var["gene_name"].tolist()
tokenized = tokenize_and_pad_batch(
    all_counts,
    np.array([vocab[g] for g in genes], dtype=int),
    max_len=len(genes) + 1,
    pad_token=vocab["<pad>"],
    pad_value=-2,
)

# 加载预训练 Transformer，分批迁移 GPU 并以 FP16 推理取 <cls> token 输出作为细胞嵌入
model = scgpt.load_pretrained(model_dir, arch_args={"n_interval": 51}).to("cuda")
model.eval()

src = torch.tensor(tokenized["genes"], dtype=torch.int64)
values = torch.tensor(tokenized["values"], dtype=torch.float32)
embs = []
with torch.no_grad(), torch.autocast("cuda", dtype=torch.float16):
    for i in range(0, src.shape[0], 64):
        s = src[i:i + 64].to("cuda")
        v = values[i:i + 64].to("cuda")
        output_dict = model(s, v, src_key_padding_mask=s.eq(vocab["<pad>"]), CLS=True)
        embs.append(output_dict["cell_embs"].float().cpu().numpy())

adata.obsm["X_scGPT"] = np.concatenate(embs, axis=0)
adata.write_h5ad("results/scgpt_integrated.h5ad")
