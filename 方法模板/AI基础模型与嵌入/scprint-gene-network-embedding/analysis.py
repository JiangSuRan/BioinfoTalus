## scPRINT 预训练单细胞大模型基因嵌入提取与细胞表征推断
import numpy as np
import scanpy as sc
import torch
from scdataloader import Collator
from scdataloader.data import SimpleAnnDataset
from torch.utils.data import DataLoader

from scprint import scPrint

model = scPrint.load_from_checkpoint("Data/scprint_checkpoint.ckpt")
adata = sc.read_h5ad("Data/cells.h5ad")

model.eval()
model.pred_log_adata = False

if "organism_ontology_term_id" not in adata.obs:
    adata.obs["organism_ontology_term_id"] = "NCBITaxon:9606"

# Collator 按模型词表组织基因子集，max_len=0 不截断（how="some" 允许部分覆盖）
adataset = SimpleAnnDataset(adata, obs_to_output=["organism_ontology_term_id"])
col = Collator(
    organisms=model.organisms,
    valid_genes=model.genes,
    max_len=0,
    how="some",
    genelist=model.genes,
)
dataloader = DataLoader(adataset, collate_fn=col, batch_size=64, num_workers=8, shuffle=False)

# 前向推断取 gene_embedding 输出 (n_cells, n_genes, d_model)，基因维取均值聚合为细胞级表征
all_embeddings = []
with torch.no_grad():
    for batch in dataloader:
        output = model(
            gene_pos=batch["genes"],
            expression=batch["x"],
            req_depth=batch["depth"],
            get_gene_emb=True,
        )
        all_embeddings.append(output["gene_embedding"].cpu().numpy())

embeddings = np.concatenate(all_embeddings, axis=0)
adata.obsm["X_scprint"] = embeddings.mean(axis=1)

adata.write_h5ad("results/scprint_embeddings.h5ad")
