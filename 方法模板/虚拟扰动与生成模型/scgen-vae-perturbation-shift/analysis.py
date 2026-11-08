## scGen 变分自编码器潜在空间向量算术跨细胞类型扰动响应预测
import scanpy as sc
import scgen

train_new = sc.read_h5ad("Data/train_kang.h5ad")

scgen.SCGEN.setup_anndata(train_new, batch_key="condition", labels_key="cell_type")
model = scgen.SCGEN(train_new)
model.train(
    max_epochs=100,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
)

# 潜在空间中刺激-对照均值差向量平移到目标细胞类型后解码，得到预测扰动表达
pred, delta = model.predict(
    ctrl_key="control",
    stim_key="stimulated",
    celltype_to_predict="CD4T",
)
pred.obs["condition"] = "pred"

model.save("results/scgen_model.pt", overwrite=True)
pred.write_h5ad("results/scgen_predicted.h5ad")
