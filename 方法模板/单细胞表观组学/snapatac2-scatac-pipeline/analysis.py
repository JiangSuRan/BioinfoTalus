## SnapATAC2 大规模单细胞 ATAC-seq 片段导入、Tile 矩阵构建、谱嵌入降维与基因活性计算
import snapatac2 as snap

fragment_files = ["Data/sample1.fragments.tsv.gz", "Data/sample2.fragments.tsv.gz"]

# 片段导入即计算 TSSe/片段数指标；多样本合并为磁盘 AnnDataSet 以承载大规模数据
adata = snap.pp.import_data(
    fragment_files,
    file="results/snap_merged.h5ads",
    chrom_sizes=snap.genome.hg38,
    min_num_fragments=500,
    sorted_by_barcode=False,
)

snap.pp.filter_cells(adata, min_tsse=10, min_counts=1000)
snap.pp.add_tile_matrix(adata, bin_size=5000)

# 黑名单区域剔除后再选 top 特征，避免 ENCODE 黑名单 bin 占据特征预算
snap.pp.select_features(adata, n_features=50000, blacklist="Data/blacklist.bed")

snap.tl.spectral(adata)
snap.pp.harmony(adata, batch="sample", key_added="X_spectral_harmony")
snap.pp.knn(adata, use_rep="X_spectral_harmony")
snap.tl.leiden(adata, resolution=0.5)
snap.tl.umap(adata, use_rep="X_spectral_harmony", use_dims=30, random_state=42)

# 基因活性矩阵：以 TSS 附近片段投入量近似基因表达
gene_matrix = snap.pp.make_gene_matrix(adata, snap.genome.hg38)
gene_matrix.write_h5ad("results/gene_activity_matrix.h5ad")
