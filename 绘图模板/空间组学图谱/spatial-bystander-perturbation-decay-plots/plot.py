## 空间微环境供体敲除分区、旁观者细胞邻近差异表达、距离衰减及特异性散点图
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

zones = pd.read_csv("Data/celcomen_bystander_zones.csv")
de_df = pd.read_csv("Data/celcomen_bystander_de.csv")

# 图1：供体敲除点位 vs 半径内/外受体细胞空间分区图
fig_z, ax_z = plt.subplots(figsize=(7.2, 6.6), constrained_layout=True)
palette = {"other": "#dddddd", "B far": "#1f77b4", "B near A (<30um)": "#ff7f0e", "A perturbed": "#d62728"}
for lbl in ["other", "B far", "B near A (<30um)", "A perturbed"]:
    sub = zones[zones["zone"] == lbl]
    ax_z.scatter(sub["x"], sub["y"], s=24 if lbl == "other" else 44, color=palette[lbl], alpha=0.88, edgecolor="none", label=lbl)
ax_z.set(xlabel="x (um)", ylabel="y (um)", title="Donor KO (A) vs spatial readout zone (B)")
ax_z.legend(loc="upper right", fontsize=8.5, framealpha=0.92)
fig_z.savefig("results/spatial_ab_zones.png", dpi=300, bbox_inches="tight")
plt.close(fig_z)

# 图2：邻域 B 细胞相对远端 B 细胞上下调 Top 差异基因双向水平条形图
de_sorted = de_df.sort_values("de_near_minus_far")
top_de = pd.concat([de_sorted.head(12), de_sorted.tail(12)])
colors = ["#3b6fb6" if v < 0 else "#c0392b" for v in top_de["de_near_minus_far"]]

fig_b, ax_b = plt.subplots(figsize=(6.8, 7.6), constrained_layout=True)
ax_b.barh(range(len(top_de)), top_de["de_near_minus_far"], color=colors, edgecolor="none")
ax_b.set_yticks(range(len(top_de)))
ax_b.set_yticklabels(top_de["gene"], fontsize=8.5)
ax_b.axvline(0, color="gray", lw=1)
ax_b.invert_yaxis()
ax_b.set(xlabel="Near-A minus far-A B-cell response", title="Bystander DE genes within KO radius")
fig_b.savefig("results/bystander_near_vs_far_de.png", dpi=300, bbox_inches="tight")
plt.close(fig_b)

# 图3：旁观者 B 细胞反事实响应随到被敲供体细胞距离的环带衰减曲线与细胞数分布
b_cells = zones[zones["is_B"]].copy()
edges = [0.0, 15.0, 30.0, 50.0, 1e9]
ring_labels = ["0-15", "15-30", "30-50", ">50"]
mags, counts = [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    sub = b_cells[(b_cells["dist_to_A"] > lo) & (b_cells["dist_to_A"] <= hi)]
    counts.append(len(sub))
    mags.append(float(sub["mean_abs_diff"].mean()) if len(sub) else np.nan)

fig_d, axes_d = plt.subplots(1, 2, figsize=(11.8, 4.6), constrained_layout=True)
axes_d[0].plot(ring_labels, mags, "o-", color="crimson", lw=2.2, ms=8.5)
axes_d[0].axvspan(-0.5, 1.5, color="#ffe0b2", alpha=0.5, label="within 30um")
axes_d[0].set(ylabel="Mean |expr change| of B cells", xlabel="Distance of B cell to nearest perturbed A (um)", title="B-cell perturbation effect vs distance to A")
axes_d[0].legend(frameon=False, fontsize=9)
axes_d[0].grid(axis="y", alpha=0.3)

axes_d[1].bar(ring_labels, counts, color="#8fb0d9", edgecolor="none")
for xi, c in enumerate(counts):
    axes_d[1].text(xi, c, str(c), ha="center", va="bottom", fontsize=9)
axes_d[1].set(ylabel="Number of B cells", xlabel="Distance ring (um)", title="B-cell counts per spatial ring")
fig_d.savefig("results/bystander_distance_decay.png", dpi=300, bbox_inches="tight")
plt.close(fig_d)

# 图4：B 细胞类型特异性旁观者 DE 与无差别全部近邻 DE 相关性对比
x = de_sorted["de_generic"].to_numpy()
y = de_sorted["de_near_minus_far"].to_numpy()
r = float(np.corrcoef(x, y)[0, 1])
lim = float(np.nanmax(np.abs(np.concatenate([x, y])))) * 1.1

fig_g, ax_g = plt.subplots(figsize=(6.2, 6.0), constrained_layout=True)
ax_g.scatter(x, y, s=16, alpha=0.55, color="#5b8c5a", edgecolor="none")
ax_g.plot([-lim, lim], [-lim, lim], "--", color="gray", lw=1)
ax_g.axhline(0, color="gray", lw=0.8)
ax_g.axvline(0, color="gray", lw=0.8)
for _, row in pd.concat([de_sorted.head(4), de_sorted.tail(4)]).iterrows():
    ax_g.annotate(row["gene"], (row["de_generic"], row["de_near_minus_far"]), fontsize=7.5)
ax_g.set(xlabel="Generic neighbor DE (all cell types)", ylabel="Cell-type B-specific DE", title=f"Cell-type-specific vs generic readout (r={r:.2f})")
fig_g.savefig("results/bystander_specific_vs_generic.png", dpi=300, bbox_inches="tight")
plt.close(fig_g)
