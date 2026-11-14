## scGPT 跨物种小鼠同源基因投影与零样本虚拟基因敲除细胞影响排序
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scgpt.model import TransformerModel
from scgpt.preprocess import binning
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import load_pretrained, set_seed
import torch

set_seed(0)
device = torch.device("cuda")

mdir = Path("Data/scgpt_human")
vocab = GeneVocab.from_file(mdir / "vocab.json")
for sp_tok in ("<pad>", "<cls>", "<eoc>"):
    if sp_tok not in vocab:
        vocab.append_token(sp_tok)
vocab.set_default_index(vocab["<pad>"])
cfg = json.load(open(mdir / "args.json"))

model = TransformerModel(
    len(vocab),
    cfg["embsize"],
    cfg["nheads"],
    cfg["d_hid"],
    cfg["nlayers"],
    nlayers_cls=cfg["n_layers_cls"],
    n_cls=1,
    vocab=vocab,
    dropout=cfg["dropout"],
    pad_token=cfg["pad_token"],
    pad_value=cfg["pad_value"],
    do_mvc=True,
    do_dab=False,
    use_batch_labels=False,
    domain_spec_batchnorm=False,
    explicit_zero_prob=False,
    use_fast_transformer=True,
    fast_transformer_backend="flash",
    pre_norm=False,
)
sd = torch.load(mdir / "best_model.pt", map_location=device)
load_pretrained(model, sd, verbose=False)
model.to(device).eval()

# 小鼠基因映射为人类同源基因（优先同源表，辅以全大写符号匹配）并对齐预训练词表
ortho_df = pd.read_csv("Data/mouse_to_human_orthologs.tsv", sep="\t").dropna()
ortho = dict(zip(ortho_df.iloc[:, 0].astype(str), ortho_df.iloc[:, 1].astype(str)))
vocab_set = set(vocab.get_stoi().keys())

adata = ad.read_h5ad("Data/mouse_cells.h5ad")
mouse_genes = adata.var_names.astype(str).tolist()
human_syms = [
    ortho[g] if (g in ortho and ortho[g] in vocab_set)
    else (g.upper() if g.upper() in vocab_set else None)
    for g in mouse_genes
]
mapped_mask = np.array([h is not None for h in human_syms])
mouse_to_human = {g: h for g, h in zip(mouse_genes, human_syms) if h is not None}

adata_h = adata[:, mapped_mask].copy()
adata_h.var_names = [h for h in human_syms if h is not None]
adata_h = adata_h[:, ~adata_h.var_names.duplicated()].copy()

# 预先对非零基因分箱并缓存序列，保证基线态与敲除态仅差目标基因移除本身
X = adata_h.X.toarray() if sp.issparse(adata_h.X) else np.asarray(adata_h.X)
X = np.asarray(X, dtype=np.float32)
gene_ids = np.array([vocab[g] for g in adata_h.var_names], dtype=np.int64)
n_bins = cfg.get("n_bins", 51)
tokens = []
for row in X:
    nz = np.nonzero(row)[0]
    if len(nz) > 1199:
        nz = np.sort(nz[np.argsort(-row[nz])[:1199]])
    tokens.append((gene_ids[nz], np.asarray(binning(row[nz], n_bins), dtype=np.float32)))

pad_id, cls_id, pad_val = vocab["<pad>"], vocab["<cls>"], cfg["pad_value"]
ko_genes = [g for g in ["Gata1", "Mpo", "Elane"] if g in mouse_to_human]

def encode_batch(drop_id):
    embs = []
    for s in range(0, len(tokens), 128):
        rows = []
        for g, v in tokens[s:s + 128]:
            if drop_id is not None:
                keep = g != drop_id
                g, v = g[keep], v[keep]
            rows.append((g, v))
        L = max(len(g) + 1 for g, _ in rows)
        G = np.full((len(rows), L), pad_id, dtype=np.int64)
        V = np.full((len(rows), L), pad_val, dtype=np.float32)
        for j, (g, v) in enumerate(rows):
            G[j, 0], V[j, 0] = cls_id, pad_val
            G[j, 1:1 + len(g)], V[j, 1:1 + len(v)] = g, v
        gt = torch.from_numpy(G).to(device)
        vt = torch.from_numpy(V).to(device)
        with torch.no_grad(), torch.autocast("cuda", dtype=torch.float16):
            e = model._encode(gt, vt, src_key_padding_mask=gt.eq(pad_id))
        embs.append(e[:, 0, :].float().cpu().numpy())
    E = np.concatenate(embs, axis=0)
    return E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-8)

base_emb = encode_batch(None)
adata_h.obsm["X_scGPT"] = base_emb
sc.pp.neighbors(adata_h, use_rep="X_scGPT", n_neighbors=15, random_state=0)
sc.tl.umap(adata_h, random_state=0)

# 单基因零样本虚拟敲除（删除目标基因 token）并按细胞类型汇总偏移量与基因总排序
rank_rows = []
first_shift = None
for mg in ko_genes:
    ko_emb = encode_batch(vocab[mouse_to_human[mg]])
    shift = np.linalg.norm(base_emb - ko_emb, axis=1)
    if first_shift is None:
        first_shift = shift
        adata_h.obs["ko_shift"] = shift
    rank_rows.append({"gene": mg, "mean_ko_impact": float(shift.mean())})

pd.DataFrame({
    "cell_id": adata_h.obs_names.astype(str),
    "cell_type": adata_h.obs["cell_type"].astype(str).values,
    "shift": first_shift,
}).to_csv("results/scgpt_mouse_ko_cell_shifts.csv", index=False)

pd.DataFrame(rank_rows).sort_values("mean_ko_impact", ascending=False).to_csv(
    "results/scgpt_mouse_ko_ranking.tsv", sep="\t", index=False
)
adata_h.write_h5ad("results/scgpt_mouse_ko_umap.h5ad")
