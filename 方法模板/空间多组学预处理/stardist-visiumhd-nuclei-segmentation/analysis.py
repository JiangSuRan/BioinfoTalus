## StarDist H&E 细胞核深度学习分割与 Visium HD 亚细胞转录本多边形归箱聚合
import geopandas as gpd
import pandas as pd
import scanpy as sc
from csbdeep.utils import normalize
from shapely.geometry import Polygon
from stardist.models import StarDist2D
import tifffile

he_img = tifffile.imread("Data/tissue_hires_image.tiff")
transcripts = pd.read_parquet("Data/detected_transcripts.parquet")

# 归一化 H&E 组织切片图像并用 2D_versatile_he 预训练权重分块预测每个细胞核多边形轮廓
img_norm = normalize(he_img, 5.0, 95.0, axis=(0, 1))
model = StarDist2D.from_pretrained("2D_versatile_he")
labels, polys = model.predict_instances_big(
    img_norm,
    axes="YXC",
    block_size=4096,
    min_overlap=128,
    context=128,
    prob_thresh=0.35,
    nms_thresh=0.3,
)

# 将多边形顶点构造为 GeoDataFrame 并通过空间包含连接（within）把单个转录本分子分配到对应细胞核
geoms = [Polygon(zip(coords[1], coords[0])) for coords in polys["coord"]]
gdf_nuclei = gpd.GeoDataFrame({"cell_id": [f"cell_{i}" for i in range(len(geoms))]}, geometry=geoms)

gdf_tx = gpd.GeoDataFrame(
    transcripts,
    geometry=gpd.points_from_xy(transcripts["x_pixel"], transcripts["y_pixel"]),
)
joined = gpd.sjoin(gdf_tx, gdf_nuclei, how="inner", predicate="within")
counts_df = joined.groupby(["cell_id", "feature_name"]).size().unstack(fill_value=0)

adata_seg = sc.AnnData(X=counts_df.to_numpy(), obs=pd.DataFrame(index=counts_df.index), var=pd.DataFrame(index=counts_df.columns))
gdf_nuclei.to_file("results/stardist_nuclei_polygons.geojson", driver="GeoJSON")
adata_seg.write_h5ad("results/stardist_segmented_cells.h5ad", compression="gzip")
