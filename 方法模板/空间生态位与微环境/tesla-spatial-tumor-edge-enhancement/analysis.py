## TESLA 组织学引导的超分辨率基因表达增强与目标细胞群轮廓标注
import cv2 as cv
import numpy as np
import scanpy as sc
import TESLA as tesla
from imageio import imread
from scipy.sparse import issparse

# 读入 Visium 计数与 HE 图像
adata = sc.read_h5ad("Data/visium_sample.h5ad")
img = imread("Data/spatial/tissue_hires_image.png")
if adata.obs["in_tissue"].dtype != bool and "x1" in adata.obs.columns:
    adata = adata[adata.obs["x1"] == 1]
adata.var_names = [i.upper() for i in list(adata.var_names)]
adata.var["genename"] = adata.var.index.astype("str")

resize_factor = 2000 / np.min(img.shape[0:2])
counts = adata
counts.var_names_make_unique()
counts.raw = counts
sc.pp.log1p(counts)
if issparse(counts.X):
    counts.X = counts.X.A

# 组织轮廓检测（Canny）与 TESLA 图像引导插值增强
cnt = tesla.cv2_detect_contour(img, apertureSize=5, L2gradient=True)
binary = np.zeros(img.shape[0:2], dtype=np.uint8)
cv.drawContours(binary, [cnt], -1, (1), thickness=-1)

enhanced_exp_adata = tesla.imputation(
    img=img,
    raw=counts,
    cnt=cnt,
    genes=counts.var.index.tolist(),
    shape="None",
    res=1,
    s=1,
    k=2,
    num_nbs=10,
)

# 目标群（T 细胞标记）在增强表达上做小目标聚类标注
TGenes = ["CD3D", "CD3E", "CD3G"]
genes = list(set([i for i in TGenes if i in enhanced_exp_adata.var.index]))
pred_refined, target_clusters, c_m = tesla.annotation(
    img=img,
    binary=binary,
    sudo_adata=enhanced_exp_adata,
    genes=genes,
    resize_factor=resize_factor,
    num_required=2,
    target_size="small",
)
topClasses = [i for i, v in c_m if v == max([v2 for _, v2 in c_m])]

ret_img = tesla.visualize_annotation(
    img=img,
    binary=binary,
    resize_factor=resize_factor,
    pred_refined=pred_refined,
    target_clusters=topClasses,
    c_m=c_m,
)
cv.imwrite("results/tesla_CD3_annotation.jpg", ret_img)

enhanced_exp_adata.write_h5ad("results/tesla_enhanced_expression.h5ad")
