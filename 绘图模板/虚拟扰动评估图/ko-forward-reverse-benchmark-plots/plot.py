## 虚拟敲除正向精度、反向识别相似度热图及细胞类型脆弱性排序图
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fwd = pd.read_csv("Data/scgpt_forward_metrics.tsv", sep="\t").sort_values("Pearson_delta", ascending=False)
delta_df = pd.read_csv("Data/scgpt_best_ko_delta.csv")
sim_df = pd.read_csv("Data/scgpt_reverse_similarity.csv", index_col=0)
shifts_df = pd.read_csv("Data/scgpt_mouse_ko_cell_shifts.csv")
rank_df = pd.read_csv("Data/scgpt_mouse_ko_ranking.tsv", sep="\t").sort_values("mean_ko_impact", ascending=False)

# 图1：正向敲除条件 Pearson(delta) 柱状图 + 最优扰动 Top DE 变化量散点
fig_f, axes_f = plt.subplots(1, 2, figsize=(12.2, 4.8), constrained_layout=True)
axes_f[0].bar(
    range(len(fwd)),
    fwd["Pearson_delta"].values,
    color=["#2ca25f" if v > 0 else "#de2d26" for v in fwd["Pearson_delta"]],
)
axes_f[0].axhline(0, color="k", lw=0.8)
axes_f[0].set_xticks(range(len(fwd)))
axes_f[0].set_xticklabels(fwd["perturbation"], rotation=45, ha="right", fontsize=8.5)
axes_f[0].set(ylabel="Pearson (delta, top DE)", title="Forward KO-effect accuracy per perturbation")

lim = float(np.abs(np.concatenate([delta_df["true_delta"], delta_df["pred_delta"]])).max()) * 1.1
axes_f[1].scatter(delta_df["true_delta"], delta_df["pred_delta"], s=42, color="#2b8cbe")
axes_f[1].plot([-lim, lim], [-lim, lim], "k--", lw=1)
axes_f[1].axhline(0, color="grey", lw=0.6)
axes_f[1].axvline(0, color="grey", lw=0.6)
axes_f[1].set_xlim(-lim, lim)
axes_f[1].set_ylim(-lim, lim)
axes_f[1].set(xlabel="True expression change", ylabel="Predicted expression change", title=f"Best KO target: {fwd.iloc[0]['perturbation']}")
fig_f.savefig("results/forward_ko_accuracy.png", dpi=300, bbox_inches="tight")
plt.close(fig_f)

# 图2：反向敲除识别相似度热图（亮对角线反映从观测表达谱精准反推真实敲除基因）
S = sim_df.to_numpy(dtype=float)
n = S.shape[0]
fig_s, ax_s = plt.subplots(figsize=(0.85 * n + 3.2, 0.78 * n + 2.5), constrained_layout=True)
im = ax_s.imshow(S, cmap="viridis", aspect="auto")
ax_s.set_xticks(range(n))
ax_s.set_yticks(range(n))
ax_s.set_xticklabels(sim_df.columns, rotation=45, ha="right", fontsize=8.5)
ax_s.set_yticklabels(sim_df.index, fontsize=8.5)
ax_s.set(xlabel="Candidate predicted KO profile", ylabel="Observed KO profile", title="Reverse inference: profile similarity")
thresh = S.max() * 0.7
for i in range(n):
    for j in range(n):
        ax_s.text(j, i, f"{S[i, j]:.2f}", ha="center", va="center", color=("white" if S[i, j] < thresh else "black"), fontsize=7.5)
fig_s.colorbar(im, ax=ax_s, shrink=0.82, label="Pearson (delta)")
fig_s.savefig("results/reverse_ko_similarity_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig_s)

# 图3：细胞类型虚拟敲除嵌入偏移箱线排序 + 多候选基因总体影响对比
cnt = shifts_df["cell_type"].value_counts()
valid_ct = cnt[cnt >= 5].index
df_valid = shifts_df[shifts_df["cell_type"].isin(valid_ct)]
ct_order = df_valid.groupby("cell_type")["shift"].median().sort_values(ascending=False).index.tolist()

fig_m, axes_m = plt.subplots(1, 2, figsize=(12.6, 5.0), constrained_layout=True)
axes_m[0].boxplot([df_valid.loc[df_valid["cell_type"] == c, "shift"].values for c in ct_order], showfliers=False)
axes_m[0].set_xticklabels(ct_order, rotation=45, ha="right", fontsize=8.5)
axes_m[0].set(ylabel="KO impact (L2 embedding shift)", title="Virtual KO vulnerability by cell type")

axes_m[1].barh(rank_df["gene"].iloc[::-1], rank_df["mean_ko_impact"].iloc[::-1], color="#3477b8")
axes_m[1].set(xlabel="Mean KO impact (embedding shift)", title="Multi-gene virtual KO impact ranking")
fig_m.savefig("results/mouse_ko_celltype_ranking.png", dpi=300, bbox_inches="tight")
plt.close(fig_m)
