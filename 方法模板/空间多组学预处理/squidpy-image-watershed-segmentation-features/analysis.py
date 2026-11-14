## Squidpy H&E 组织影像容器构建、高斯平滑、分水岭核分割与 Spot 形态特征提取
import pandas as pd
import scanpy as sc
import squidpy as sq

adata = sc.read_h5ad("Data/visium_with_hires_image.h5ad")
lib_id = list(adata.uns["spatial"].keys())[0]
img_arr = adata.uns["spatial"][lib_id]["images"]["hires"]
scale_hires = adata.uns["spatial"][lib_id]["scalefactors"]["tissue_hires_scalef"]

# 构建 ImageContainer 并对齐 hires 缩放比例，高斯平滑后执行分水岭分割并聚合至 Spot 尺度特征
img = sq.im.ImageContainer(img_arr, layer="image", scale=scale_hires)
sq.im.process(img=img, layer="image", method="smooth", sigma=1, layer_added="smooth")
sq.im.segment(img=img, layer="smooth", method="watershed", channel=0, layer_added="segmented_watershed")

feat_summary = sq.im.calculate_image_features(
    adata,
    img,
    layer="image",
    features=["summary", "segmentation"],
    features_kwargs={"segmentation": {"label_layer": "segmented_watershed"}},
    key_added="img_features",
    n_jobs=1,
    show_progress_bar=False,
    copy=False,
)

img_feat_df = adata.obsm["img_features"]
img_feat_df.to_csv("results/squidpy_spot_image_features.csv")
adata.write_h5ad("results/squidpy_segmented_features.h5ad", compression="gzip")
