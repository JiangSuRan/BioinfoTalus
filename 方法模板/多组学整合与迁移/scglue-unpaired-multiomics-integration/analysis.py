## scGLUE 非配对单细胞转录组（scRNA-seq）与染色质可及性（scATAC-seq）先验引导图整合及标签迁移
from itertools import chain

import anndata as ad
import scglue

rna = ad.read_h5ad("Data/scrna.h5ad")
atac = ad.read_h5ad("Data/scatac.h5ad")

# GTF 注释补齐基因坐标，以基因为锚构建基因-peak 先验引导图（峰按坐标挂接到 TSS 邻近基因）
scglue.data.get_gene_annotation(rna, gtf="Data/V3_genes.gtf", gtf_by="gene_name")
guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
scglue.graph.check_graph(guidance, [rna, atac])

# ATAC 侧 LSI 降维为模型输入表征；两侧均以 NB 似然配 counts 层与批次列
scglue.data.lsi(atac, n_components=100, n_iter=15, use_highly_variable=True)

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True, use_layer="counts", use_rep="X_pca", use_batch="sample"
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True, use_rep="X_lsi", use_batch="sample"
)

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index,
)).copy()

glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    init_kws={"h_dim": 512, "random_seed": 666},
    fit_kws={"directory": "results/glue_model", "data_batch_size": 8192},
)
glue.save("results/glue_model.dill")

# 共享 GLUE 空间提取联合嵌入，并以 kNN 迁移 RNA 侧细胞类型注释到 ATAC
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
scglue.data.transfer_labels(rna, atac, "cell_type_l4", use_rep="X_glue", n_jobs=-1)

atac.write_h5ad("results/atac_glue_annotated.h5ad")
