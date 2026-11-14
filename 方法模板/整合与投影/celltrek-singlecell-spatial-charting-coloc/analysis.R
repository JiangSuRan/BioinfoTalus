## CellTrek 单细胞空间反投射精确定位（Charting）与空间共定位/邻近距离分析
library(Seurat)
library(CellTrek)
library(dplyr)

st_obj <- readRDS("Data/visium_spatial.rds")
sc_obj <- readRDS("Data/scrna_reference.rds")

# traint 构建 scRNA-seq 与 Visium 空间转录组共同嵌入并插值空间坐标后映射到单细胞分辨率
traint_obj <- traint(st_data = st_obj, sc_data = sc_obj, st_assay = "Spatial", sc_assay = "RNA", norm = "LogNormalize", nfeatures = 2000)
celltrek_res <- celltrek(
  st_sc_int = traint_obj,
  int_assay = "traint",
  sc_data = sc_obj,
  sc_assay = "RNA",
  reduction = "pca",
  intp = TRUE,
  intp_pnt = 5000,
  intp_lin = FALSE,
  nPCs = 30,
  ntree = 1000,
  dist_thresh = 0.55,
  top_spot = 5,
  spot_n = 5,
  repel_r = 20,
  repel_iter = 20,
  keep_model = FALSE
)$celltrek

# SColoc 基于 Kullback-Leibler 散度量化空间重排后细胞亚群间的显著共定位/互斥拓扑
scoloc_res <- scoloc(celltrek_res, col_cell = "celltype", use_method = "KL", eps = 1e-20)
kdist_res <- kdist(inp_df = celltrek_res@meta.data, ref = "Tumor", que = "T_cell", k = 15, new_name = "dist_to_tumor")

saveRDS(celltrek_res, "results/celltrek_mapped_seurat.rds")
saveRDS(scoloc_res, "results/celltrek_scoloc.rds")
write.csv(kdist_res, "results/celltrek_kdist.csv", row.names = TRUE)
