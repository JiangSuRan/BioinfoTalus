## STATE (Arc Institute) 单细胞状态基础模型嵌入推理
from state.emb.inference import Inference

inf = Inference()
inf.load_model("Data/state_model.ckpt")

# 批量编码：嵌入写回 obsm 并导出（键名默认 X_emb）
inf.encode_adata(
    input_adata_path="Data/cells.h5ad",
    output_adata_path="results/state_emb.h5ad",
    emb_key="X_state",
)
