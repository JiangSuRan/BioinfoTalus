## scCODA 效应热图（可信效应按 Inclusion probability 过滤）
import glob
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plot_order = {"cross_area_subclass": ["L2/3 IT", "L4 IT", "L5 IT", "L5 PT", "L6 CT", "L6 IT"]}
plot_size = {"cross_area_subclass": (-1.5, 1.5)}
covariate_order = ["sex", "region"]
drop_covariates = []

results_tables = glob.glob(os.path.join("Data", "sccoda", "**", "*_results.csv"), recursive=True)
results_table = None
for z in results_tables:
    tmp = pd.read_csv(z, index_col=0)
    tmp["Covariate"] = [
        re.sub("C\\([^,]+, Treatment\\('[^\\)]+'\\)\\)\\[T\\.([^\\]]+)\\]$", "\\1", y)
        for y in tmp["Covariate"]
    ]
    results_table = tmp if results_table is None else pd.concat([results_table, tmp], axis=0)

# Final Parameter != 0 即可信效应；Inclusion probability < 0.8 的效应置零
results_table["Credible"] = results_table["Final Parameter"] != 0
results_table["Cell Type"] = pd.Categorical(results_table["Cell Type"], categories=plot_order["cross_area_subclass"])
results_table["Covariate"] = pd.Categorical(results_table["Covariate"], categories=covariate_order)
results_table = results_table.loc[[l not in drop_covariates for l in results_table["Covariate"]], :]

tmp = results_table.loc[:, ["Covariate", "Cell Type", "Final Parameter"]].groupby(["Covariate", "Cell Type"], observed=True).mean().reset_index()
tmp.columns = ["Covariate", "Cell Type", "Effect Size"]

tmp2 = results_table.loc[:, ["Covariate", "Cell Type", "Inclusion probability"]].groupby(["Covariate", "Cell Type"], observed=True).mean().reset_index()
tmp.loc[tmp2["Inclusion probability"] < 0.8, "Effect Size"] = 0

fig, ax = plt.subplots(figsize=(8, 4))
sns.heatmap(
    tmp.pivot(index="Covariate", columns="Cell Type", values="Effect Size"),
    center=0,
    cmap="RdBu_r",
    vmin=plot_size["cross_area_subclass"][0],
    vmax=plot_size["cross_area_subclass"][1],
    xticklabels=True,
    ax=ax,
)
fig.savefig("results/sccoda_effect_size_credible.pdf", bbox_inches="tight")
tmp.to_csv("results/sccoda_effect_size_p_gt_0.8.csv")
