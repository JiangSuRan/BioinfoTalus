## scGPT TransformerGenerator 正向基因敲除预测精度评估与反向扰动基因识别
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import scipy.sparse as sp
import torch
from gears.utils import create_cell_graph_dataset_for_prediction
from scgpt.model import TransformerGenerator
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed
from torch_geometric.loader import DataLoader as GDataLoader

set_seed(42)
device = torch.device("cuda")

full = ad.read_h5ad("Data/perturb_processed.h5ad")
cond_name = dict(full.obs[["condition", "condition_name"]].values)
conds = [c for c in full.obs["condition"].unique() if c != "ctrl"][:8]

# 构建信息基因面板：强制保留全部条件的 top-20 差异基因及敲除靶基因，其余按方差补齐
force_idx = set()
name_to_var = {nm: lab for lab, nm in zip(full.var.index, full.var["gene_name"])}
for c in conds:
    for r in full.uns["top_non_dropout_de_20"][cond_name[c]]:
        force_idx.add(full.var.index.get_loc(r))
    for g in c.split("+"):
        if g != "ctrl" and g in name_to_var:
            force_idx.add(full.var.index.get_loc(name_to_var[g]))

X_full = full.X.toarray() if sp.issparse(full.X) else np.asarray(full.X)
panel = list(force_idx)
for i in np.argsort(-X_full.var(axis=0)):
    if len(panel) >= 5060:
        break
    if int(i) not in force_idx:
        panel.append(int(i))

adata = full[:, sorted(panel)].copy()

model_dir = Path("Data/scgpt_human")
vocab = GeneVocab.from_file(model_dir / "vocab.json")
for sp_tok in ("<pad>", "<cls>", "<eoc>"):
    if sp_tok not in vocab:
        vocab.append_token(sp_tok)
vocab.set_default_index(vocab["<pad>"])
cfg = json.load(open(model_dir / "args.json"))

model = TransformerGenerator(
    len(vocab),
    cfg["embsize"],
    cfg["nheads"],
    cfg["d_hid"],
    cfg["nlayers"],
    nlayers_cls=cfg["n_layers_cls"],
    n_cls=1,
    vocab=vocab,
    dropout=0,
    pad_token="<pad>",
    pad_value=0,
    pert_pad_id=0,
    use_fast_transformer=True,
)
sd = torch.load("Data/scgpt_finetuned_perturbation.pt", map_location=device)
msd = model.state_dict()
msd.update({k: v for k, v in sd.items() if k in msd and v.shape == msd[k].shape})
model.load_state_dict(msd)
model.to(device).eval()

genes = adata.var["gene_name"].tolist()
gene_ids = np.array([vocab[g] if g in vocab else vocab["<pad>"] for g in genes], dtype=int)
ctrl = adata[adata.obs["condition"] == "ctrl"]
X_panel = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
ctrl_mean = np.asarray(ctrl.X.toarray() if sp.issparse(ctrl.X) else ctrl.X).mean(axis=0)
g2i = {g: i for i, g in enumerate(adata.var.index)}

rows, pred_delta, true_delta = [], {}, {}
for c in conds:
    pert_genes = [g for g in c.split("+") if g != "ctrl"]
    cg = create_cell_graph_dataset_for_prediction(pert_genes, ctrl, genes, device, num_samples=200)
    loader = GDataLoader(cg, batch_size=64, shuffle=False)
    preds = []
    for batch in loader:
        with torch.no_grad():
            preds.append(model.pred_perturb(batch, "all", gene_ids=gene_ids, amp=True).cpu())
    pm = torch.cat(preds, dim=0).float().numpy().mean(axis=0)
    pred_delta[c] = pm - ctrl_mean

    c_mask = (adata.obs["condition"] == c).to_numpy()
    tm = X_panel[c_mask].mean(axis=0)
    true_delta[c] = tm - ctrl_mean
    de = [g2i[r] for r in full.uns["top_non_dropout_de_20"][cond_name[c]] if r in g2i][:20]
    rows.append({
        "perturbation": c,
        "Pearson_expr": pearsonr(pm[de], tm[de])[0],
        "Pearson_delta": pearsonr(pred_delta[c][de], true_delta[c][de])[0],
        "dir_acc": float(np.mean(np.sign(pred_delta[c][de]) == np.sign(true_delta[c][de]))),
        "MSE_delta": float(np.mean((pred_delta[c][de] - true_delta[c][de]) ** 2)),
    })

res = pd.DataFrame(rows).sort_values("Pearson_delta", ascending=False).reset_index(drop=True)
res.to_csv("results/scgpt_forward_metrics.tsv", sep="\t", index=False)

best = res.iloc[0]["perturbation"]
best_de = [g2i[r] for r in full.uns["top_non_dropout_de_20"][cond_name[best]] if r in g2i][:20]
pd.DataFrame({
    "gene": adata.var.index[best_de].astype(str),
    "true_delta": true_delta[best][best_de],
    "pred_delta": pred_delta[best][best_de],
}).to_csv("results/scgpt_best_ko_delta.csv", index=False)

# 反向敲除基因识别：用所有条件 Top-DE 并集构造查询谱-预测谱相似度矩阵
de_union = sorted({
    g2i[r] for c in conds for r in full.uns["top_non_dropout_de_20"][cond_name[c]] if r in g2i
})
S = np.zeros((len(conds), len(conds)))
for i, q in enumerate(conds):
    for j, cand in enumerate(conds):
        S[i, j] = pearsonr(true_delta[q][de_union], pred_delta[cand][de_union])[0]

pd.DataFrame(S, index=conds, columns=conds).to_csv("results/scgpt_reverse_similarity.csv")
