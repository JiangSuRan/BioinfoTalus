## stLearn SME 空间归一化、配受体 CCI 检验与细胞类型水平通讯网络
import matplotlib.pyplot as plt
import pandas as pd
import stlearn as st

# 读取 Visium 与 HE 图像（SME 邻域加权依赖图像形态）
data = st.Read10X("Data/visium_sample/outs")
data.var_names_make_unique()
st.add.image(adata=data, imgpath="Data/spatial/tissue_hires_image.png",
             library_id="sample", visium=True)

st.pp.filter_genes(data, min_cells=3)
st.pp.normalize_total(data)

# 上游反卷积主导标签与组成矩阵写入 obs/uns（键名一致供 CCI 使用）
spot_mixtures = pd.read_csv("Data/spotlight_deconvolution_weights.csv", index_col=0)
data.obs["cell_type"] = spot_mixtures.loc[:, "spotlight.cluster"].values.astype(str)
data.obs["cell_type"] = data.obs["cell_type"].astype("category")
data.uns["cell_type"] = spot_mixtures.drop(["spotlight.cluster"], axis=1)

# NATMI 配受体库 + spot 邻域 CCI 检验（FDR 校正）+ 细胞类型水平汇总
lrs = st.tl.cci.load_lrs(["connectomeDB2020_lit"], species="human")
st.tl.cci.run(
    data, lrs,
    min_spots=20,
    distance=None,
    n_pairs=1000,
    n_cpus=8,
)
st.tl.cci.adj_pvals(data, correct_axis="spot", pval_adj_cutoff=0.05, adj_method="fdr_bh")
st.tl.cci.run_cci(data, "cell_type", min_spots=3, spot_mixtures=True, cell_prop_cutoff=0.2, sig_spots=True, n_perms=100)

fig = st.pl.ccinet_plot(data, "cell_type", return_pos=False)
fig.savefig("results/stlearn_cci_network.png", dpi=300, bbox_inches="tight")
plt.close("all")

data.filename = "results/stlearn_cci.h5ad"
