## scPRINT 细胞类型定向基因调控网络推断与候选基因关联邻居排序
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scdataloader import Preprocessor
from scprint import scPrint
from scprint.tasks import GNInfer

adata = sc.read_h5ad("Data/retina_counts.h5ad")
ckpt_path = "Data/scprint_checkpoint.ckpt"
meta = torch.load(ckpt_path, map_location="cuda", weights_only=False)
model = scPrint.load_from_checkpoint(
    ckpt_path,
    precpt_gene_emb=None,
    classes=meta["hyper_parameters"]["label_counts"],
).cuda().eval()

# 数据经 scdataloader 统一预处理后，在目标细胞类型中取候选基因并补充高表达背景以推断局部调控网络
adata = Preprocessor(force_preprocess=True, skip_validate=True, do_postp=False)(adata)
cell_type = "retinal rod cell"
cand_symbols = ["CRX", "NRL", "NR2E3", "VSX2", "RHO", "RCVRN", "PDE6B", "GNAT1"]

sym_series = adata.var["feature_name"].astype(str)
sym_to_ens = dict(zip(sym_series, adata.var_names.astype(str)))
cand_ids = [sym_to_ens[s] for s in cand_symbols if s in sym_to_ens]

mask = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
X_sub = adata.X[mask]
mean_expr = np.asarray(X_sub.mean(axis=0)).ravel() if sp.issparse(X_sub) else np.asarray(X_sub).mean(axis=0)
bg_ids = adata.var_names.astype(str)[np.argsort(mean_expr)[::-1][:120]].tolist()
gene_set = list(dict.fromkeys(cand_ids + bg_ids))

gn = GNInfer(
    how="given",
    genes=gene_set,
    num_genes=len(gene_set),
    max_len=len(gene_set),
    max_cells=512,
    batch_size=64,
    num_workers=4,
    dtype=torch.float16,
    precision="16-mixed",
    cell_type_col="cell_type",
    head_agg="mean",
    filtration="thresh",
)
grn = gn(model, adata, cell_type=cell_type)

mat = grn.varp["GRN"]
mat = mat.toarray() if sp.issparse(mat) else np.asarray(mat)
gene_names = grn.var["feature_name"].astype(str).to_numpy()
pd.DataFrame(mat, index=gene_names, columns=gene_names).to_csv("results/scprint_grn_matrix.csv")

# 取正反两向注意力的最大绝对值作为无向关联强度，为每个候选基因提取前 15 个关联伙伴
ens_to_idx = {e: i for i, e in enumerate(grn.var_names.astype(str))}
undirected = np.maximum(np.abs(mat), np.abs(mat).T)
rows = []
for ce in cand_ids:
    if ce not in ens_to_idx:
        continue
    i = ens_to_idx[ce]
    w = undirected[i].copy()
    w[i] = 0.0
    rank = 0
    for j in np.argsort(w)[::-1]:
        if w[j] <= 0 or rank >= 15:
            break
        rank += 1
        rows.append({
            "candidate": str(gene_names[i]),
            "neighbor": str(gene_names[j]),
            "weight": float(w[j]),
            "rank": rank,
        })

pd.DataFrame(rows).to_csv("results/scprint_candidate_neighbors.csv", index=False)
