## 基于空间距离图的肿瘤/免疫细胞聚落（Colony）连通域识别
import cv2 as cv
import numpy as np
import pandas as pd
from math import ceil, pi, sqrt

cell_tab = pd.read_csv("Data/codex_celltype_area.csv", header=None)

# 图幅按坐标范围向上取整并外扩 200 px，保证边缘细胞入图
x_max = float(cell_tab[3].max())
y_max = float(cell_tab[4].max())
canvas_shape = (
    int(100 * ceil(round(y_max / 100, 1)) + 200),
    int(100 * ceil(round(x_max / 100, 1)) + 200),
)

records = []
for ct in cell_tab[6].unique():
    sub = cell_tab[cell_tab[6] == ct]
    # 细胞绘制半径由该类型面积中位数推算，使点图上细胞尺度接近真实
    r = round(sqrt(np.median(sub[5]) / pi), 1)
    image = np.zeros(canvas_shape, dtype=np.uint8)
    for _, row in sub.iterrows():
        cv.circle(image, (int(round(row[3])), int(round(row[4]))), int(r), 255, -1)

    # OTSU 二值化后以椭圆核膨胀桥接邻近细胞，物理连续的同型细胞群融合为连通域
    ret, thresh = cv.threshold(image, 9, 255, 0, cv.THRESH_OTSU)
    kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (30, 30))
    dilated = cv.dilate(thresh, kernel, iterations=1)
    contours, hierarchy = cv.findContours(dilated, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)

    # 每个细胞归入其质心所在的连通域，连通域索引即聚落编号
    for _, row in sub.iterrows():
        pt = (int(round(row[3])), int(round(row[4])))
        colony_id = -1
        for k, contour in enumerate(contours):
            if cv.pointPolygonTest(contour, pt, False) >= 0:
                colony_id = k
                break
        records.append({"x": row[3], "y": row[4], "area": row[5], "celltype": row[6], "colony_id": colony_id})

labels = pd.DataFrame(records)
labels.to_csv("results/spatial_colony_labels.csv", index=False)
