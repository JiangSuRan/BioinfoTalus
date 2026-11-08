## topact 亚细胞空间转录组密度场构建与 blob 拓扑峰检测
import numpy as np
import pandas as pd
from topact.densetools import density

# Stereo-seq bin20 细胞分类表（坐标、样本、簇标签）
bin20 = pd.read_csv(
    "Data/bin20_classifications.csv",
    usecols=["x", "y", "mouse_ID", "clusterid_0_6"],
)

# 取目标簇（肾小球足细胞 pods）的空间分布
bin20pods = bin20[bin20["clusterid_0_6"].str.match(r"(.)*pods")]

offsets_csv = pd.read_csv("Data/offsets.csv")
offsets = {}
for row in offsets_csv.itertuples():
    offsets[row.sample] = (row.xmin, row.ymin)

# 按样本构建 20 像素分格指示矩阵，密度场平滑后 blob_dog 检测肾小球团块中心
for sample in ["control_1", "control_2", "4032_1", "4030_1"]:
    df = bin20pods[bin20pods.mouse_ID == sample]
    xmin, xmax, ymin, ymax = df.x.min(), df.x.max(), df.y.min(), df.y.max()
    N, M = xmax - xmin + 1, ymax - ymin + 1

    matrix = np.zeros((M // 20 + 1, N // 20 + 1))
    for row in df.itertuples():
        matrix[(row.y - ymin) // 20, (row.x - xmin) // 20] = 1

    d = density(matrix, 1)

    pods = skimage_blob_dog(d)

    col = pods[:, 0].copy()
    pods[:, 0] = pods[:, 1] + xmin
    pods[:, 1] = col + ymin

    np.savetxt(f"results/{sample}_gloms.txt", pods, fmt="%1.f")


def skimage_blob_dog(d):
    import skimage
    blobs = skimage.feature.blob_dog(d, threshold=0.2)[:, :2] * 20 + 10
    return blobs
