## Drug2Cell 基于 ChEMBL 药靶互作数据库的单细胞药物靶点富集与敏感性打分
import drug2cell as d2c
import scanpy as sc

adata = sc.read_h5ad("Data/scrna_preprocessed.h5ad")

# 按 ChEMBL 药物-靶基因映射计算单细胞药物活动评分（存放于 adata.uns['drug2cell'] 兼容对象中）
d2c.score(adata, use_raw=False)

# 对高评分药物执行亚群 Wilcoxon 检验及药靶超几何富集，识别特定细胞群敏感药物分子
sc.tl.rank_genes_groups(adata.uns["drug2cell"], method="wilcoxon", groupby="cell_type")
overrepresented = d2c.hypergeometric(adata)

drug_diff = sc.get.rank_genes_groups_df(adata.uns["drug2cell"], group=None)
drug_diff.to_csv("results/drug2cell_cluster_drug_scores.csv", index=False)
overrepresented.to_csv("results/drug2cell_hypergeometric_enrichment.csv")
adata.write_h5ad("results/drug2cell_scored.h5ad", compression="gzip")
