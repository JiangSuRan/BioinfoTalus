## Scirpy 单细胞免疫组库 VDJ 用途、CDR3 谱系长度、样本重叠热图与克隆型网络
import matplotlib.pyplot as plt
import scanpy as sc
import scirpy as ir

adata = sc.read_h5ad("Data/scirpy_tcr_annotated.h5ad")

# VDJ 基因片段使用堆叠图、CDR3 氨基酸长度频谱图与样本间克隆重叠系数矩阵联合刻画组库偏倚
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
ir.pl.vdj_usage(adata, full_combination=False, top_n=12, ax=axes[0, 0])
ir.pl.spectratype(adata, color="cell_type", viztype="curve", ax=axes[0, 1])
ir.pl.repertoire_overlap(adata, groupby="sample_id", ax=axes[1, 0])
ir.pl.clonotype_network(adata, color="cell_type", base_size=10, label_fontsize=8, ax=axes[1, 1])
fig.tight_layout()
fig.savefig("results/scirpy_repertoire_visual_suite.pdf", dpi=300, bbox_inches="tight")
