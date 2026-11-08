## PASTE 融合 Gromov-Wasserstein 最优传输空间切片两两对齐与坐标重排
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
import paste as pst
import squidpy as sq

def load_and_preprocess(sample_path):
    visium_sample = sq.read.visium(sample_path)
    sc.pp.filter_genes(visium_sample, min_counts=15)
    sc.pp.filter_cells(visium_sample, min_counts=100)
    visium_sample.var_names_make_unique()
    return visium_sample

slice1 = load_and_preprocess("Data/slice1/outs")
slice2 = load_and_preprocess("Data/slice2/outs")
slices = [slice1, slice2]

# alpha 加权表达 KL 与空间 GW 成本；PASTE2 局部对齐返回软匹配矩阵 pi
alpha = 0.1
cost = "kl"

slice1.obsm["rgb"] = slice1.obsm["spatial"]
slice2.obsm["rgb"] = slice2.obsm["spatial"]

from paste2.PASTE2 import partial_pairwise_align

slice12_al = partial_pairwise_align(slice1, slice2, s=alpha, dissimilarity=cost)
np.save("results/paste2_partial_pairwise.npy", slice12_al)

# 按 pi 堆叠两切片得到对齐后坐标
new_slices = pst.stack_slices_pairwise(slices, [slice12_al])

slice_colors = ["#e41a1c", "#377eb8"]
fig = plt.figure(figsize=(7, 7))
for i in range(len(new_slices)):
    pst.plot_slice(new_slices[i], slice_colors[i], s=70)
plt.legend(handles=[
    mpatches.Patch(color=slice_colors[0], label="slice1"),
    mpatches.Patch(color=slice_colors[1], label="slice2"),
])
plt.gca().invert_yaxis()
plt.axis("off")
fig.savefig("results/paste2_alignment.pdf", dpi=300, bbox_inches="tight")

# 导出对齐前后坐标表（与 tissue_positions 列格式一致）
sample_name_list = ["slice1", "slice2"]
for i in range(2):
    df_tissue = slices[i].obs[["in_tissue", "array_row", "array_col"]]
    df_coord_orig = pd.DataFrame(slices[i].obsm["spatial"], index=df_tissue.index,
                                 columns=["pxl_row_in_fullres", "pxl_col_in_fullres"])
    df_coord_new = pd.DataFrame(new_slices[i].obsm["spatial"], index=df_tissue.index,
                                columns=["pxl_row_in_fullres_new", "pxl_col_in_fullres_new"])
    df_coord = pd.concat([df_tissue, df_coord_orig, df_coord_new], axis=1).reset_index().rename(columns={"index": "barcode"})
    df_coord.to_csv(f"results/paste2_new_tissue_positions_{sample_name_list[i]}.csv", index=False)
