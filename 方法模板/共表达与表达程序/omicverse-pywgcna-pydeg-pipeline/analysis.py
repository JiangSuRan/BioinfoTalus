## omicverse bulk 差异表达（pyDEG）与加权共表达模块（pyWGCNA）
import numpy as np
import pandas as pd
import omicverse as ov

counts = pd.read_csv("Data/bulk_counts.csv", index_col=0)

# pyDEG：DESeq2 尺寸因子归一化差异表达
deg = ov.bulk.pyDEG(counts)
deg.drop_duplicates_index()
deg_result = deg.deseq2(
    test_list=counts.columns[:3].tolist(),
    control_list=counts.columns[3:6].tolist(),
)
deg_result = deg.foldchange_set(fc_threshold=-1, pval_threshold=0.05, logp_max=6)

deg_result.to_csv("results/ov_pydeg_results.csv")

# pyWGCNA：样本 × 基因表达矩阵构建共表达模块
data = counts.T
wgcna = ov.bulk.pyWGCNA(name="bulk_analysis", species="human", geneExp=data, save=True)
wgcna.preprocess()
wgcna.findModules()
wgcna.analyseWGCNA()

wgcna.save("results/ov_pywgcna")
module_genes = wgcna.get_module_genes()
module_genes.to_csv("results/ov_pywgcna_module_genes.csv")
