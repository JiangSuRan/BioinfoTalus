## scvi-tools (SCVI) 深度生成模型单细胞多批次概率整合与潜在表征提取
import scanpy as sc
import scvi

adata = sc.read_h5ad("Data/counts_batched.h5ad")

# seurat_v3 HVG 基于计数层且按批次分组估计，避免批次间深度差主导特征选择
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="batch",
    subset=True,
)

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata, n_latent=20, gene_likelihood="nb", encode_covariates=True)
model.train(max_epochs=500)

adata.obsm["X_scVI"] = model.get_latent_representation()

model.save("results/scvi_model/", overwrite=True)
adata.write_h5ad("results/scvi_integrated.h5ad")
