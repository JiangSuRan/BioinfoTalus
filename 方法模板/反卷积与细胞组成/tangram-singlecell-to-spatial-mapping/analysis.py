## Tangram 空间转录组单细胞映射、基因表达空间填补与细胞类型概率投影
import scanpy as sc
import tangram as tg

adata_sc = sc.read_h5ad("Data/scrna_annotated.h5ad")
adata_sp = sc.read_h5ad("Data/visium_spatial.h5ad")

sc.tl.rank_genes_groups(adata_sc, groupby="cell_type", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata_sc, group=None).groupby("group").head(40)["names"].unique().tolist()

# 筛选 scRNA-seq 与空间数据共有训练标志基因并统一排序索引，构造非负余弦相似度优化目标
tg.pp_adatas(adata_sc, adata_sp, genes=markers)
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_sp,
    mode="cells",
    density_prior="rna_count_based",
    num_epochs=500,
    device="cpu",
)

tg.project_cell_annotations(ad_map, adata_sp, annotation="cell_type")
ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)

ad_map.write_h5ad("results/tangram_cell_to_spot_map.h5ad", compression="gzip")
adata_sp.write_h5ad("results/tangram_spatial_annotated.h5ad", compression="gzip")
ad_ge.write_h5ad("results/tangram_imputed_genes.h5ad", compression="gzip")
