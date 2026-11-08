## Nicheformer 空间与解离单细胞基础模型零样本特征嵌入提取
import numpy as np
import scanpy as sc
import torch

from nicheformer.data.dataset import tokenize_data
from nicheformer.models._nicheformer import Nicheformer

adata = sc.read_h5ad("Data/spatial_dissociated.h5ad")
model = Nicheformer.load_from_checkpoint("Data/nicheformer_checkpoint.ckpt")
model.eval()

# 基因表达按每基因中位计数缩放后 token 化（4096 序列位 + 辅助 token）
X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
median_counts = np.array([np.median(col[col > 0]) if np.any(col > 0) else 1.0 for col in X.T])
tokens = tokenize_data(X, median_counts_per_gene=median_counts).astype(np.int64)

# 分批前向，transformer 输出均值池化为 512 维细胞嵌入
embeddings = []
with torch.no_grad():
    for i in range(0, tokens.shape[0], 256):
        x = torch.tensor(tokens[i:i + 256])
        attention_mask = torch.ones_like(x, dtype=torch.bool)
        output = model(x, attention_mask)
        embeddings.append(output["transformer_output"].mean(dim=1).cpu().numpy())

adata.obsm["X_nicheformer"] = np.concatenate(embeddings, axis=0)
adata.write_h5ad("results/nicheformer_emb.h5ad")
