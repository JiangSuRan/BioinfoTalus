## 生成式扰动 PPI 去偏差置信区间、火山图、邻居偏差一致性与真实预测相关性评估图
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr

df = pd.read_csv("Data/ppi_gene_summary.csv")
bias_df = pd.read_csv("Data/ppi_neighbor_bias.csv", index_col=0)

# 图1：Top 效应基因观测 vs PPI 去偏差 95% 置信区间横向误差棒对比
top_df = df.sort_values("abs_effect", ascending=False).head(20).reset_index(drop=True)
z_crit = norm.ppf(0.975)
y = np.arange(len(top_df))

fig_ci, ax_ci = plt.subplots(figsize=(10.5, 8.0), constrained_layout=True)
ax_ci.barh(
    y + 0.2, top_df["obs_mean"], height=0.36, xerr=z_crit * top_df["obs_se"],
    color="#1f77b4", alpha=0.85, label="Observed (real)",
    error_kw=dict(ecolor="#16324f", lw=1.0, capsize=3),
)
ax_ci.barh(
    y - 0.2, top_df["ppi_mean"], height=0.36, xerr=z_crit * top_df["ppi_se"],
    color="#ff7f0e", alpha=0.85, label="Debiased (PPI)",
    error_kw=dict(ecolor="#8a4500", lw=1.0, capsize=3),
)
ax_ci.set_yticks(y)
ax_ci.set_yticklabels(top_df["gene"].astype(str), fontsize=8.5)
ax_ci.set_xlabel("Mean gene expression (log-normalized)")
ax_ci.set_title("Top 20 effect-size genes: 95% CI (observed vs debiased)")
ax_ci.legend(loc="lower right", fontsize=9.5)
ax_ci.invert_yaxis()
ax_ci.margins(y=0.01)
fig_ci.savefig("results/ppi_confidence_intervals.png", dpi=300, bbox_inches="tight")
plt.close(fig_ci)

# 图2：PPI 去偏差差异表达火山图（Bonferroni 阈值）
ppi_eff = (df["ppi_mean"] - df["ctrl_mean"]).to_numpy()
obs_eff = (df["obs_mean"] - df["ctrl_mean"]).to_numpy()
neg_log_p = -np.log10(df["pval"].to_numpy() + 1e-300)
sig_thr = 0.05 / len(df)
sig = df["pval"].to_numpy() < sig_thr
up = sig & (obs_eff > 0)
dn = sig & (obs_eff < 0)

fig_v, ax_v = plt.subplots(figsize=(8.2, 6.4), constrained_layout=True)
ax_v.scatter(ppi_eff[~sig], neg_log_p[~sig], c="lightgrey", alpha=0.45, s=12, label="Not significant", edgecolor="none")
ax_v.scatter(ppi_eff[up], neg_log_p[up], c="#d62728", alpha=0.78, s=22, label=f"Sig. & up (n={int(up.sum())})", edgecolor="none")
ax_v.scatter(ppi_eff[dn], neg_log_p[dn], c="#1f77b4", alpha=0.78, s=22, label=f"Sig. & down (n={int(dn.sum())})", edgecolor="none")
ax_v.axhline(-np.log10(sig_thr), color="grey", linestyle="--", lw=1.0, alpha=0.7)
ax_v.axvline(0.0, color="grey", linestyle=":", lw=0.8, alpha=0.6)
ax_v.set(xlabel="Debiased effect size (PPI mean - control)", ylabel="-log10(p-value)", title="Differential expression volcano (PPI debiased)")
ax_v.legend(loc="upper center", fontsize=9)
fig_v.savefig("results/ppi_debiased_volcano.png", dpi=300, bbox_inches="tight")
plt.close(fig_v)

# 图3：同簇参考扰动偏差向量一致性散点（验证系统偏差跨扰动可迁移性）
nbs = bias_df.columns.tolist()
base = nbs[0]
others = nbs[1:min(4, len(nbs))]
fig_b, axes_b = plt.subplots(1, len(others), figsize=(4.5 * len(others), 4.5), constrained_layout=True, squeeze=False)
axes_flat = axes_b.ravel()

for i, other in enumerate(others):
    ax = axes_flat[i]
    b1 = bias_df[base].to_numpy()
    b2 = bias_df[other].to_numpy()
    r = pearsonr(b1, b2)[0]
    lim = float(max(np.abs(b1).max(), np.abs(b2).max())) * 1.05
    ax.scatter(b1, b2, alpha=0.4, s=14, c="steelblue", edgecolor="none")
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=1.0, alpha=0.7)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set(xlabel=f"{base} bias", ylabel=f"{other} bias", title=f"Within-cluster bias (r = {r:.3f})")

fig_b.suptitle("Bias-consistency check across neighbor perturbations", fontsize=12)
fig_b.savefig("results/ppi_bias_consistency.png", dpi=300, bbox_inches="tight")
plt.close(fig_b)

# 图4：真实 vs 生成平均表达与效应量双联相关性散点
obs_mean = df["obs_mean"].to_numpy()
gen_mean = df["ppi_mean"].to_numpy()
r_mean = pearsonr(obs_mean, gen_mean)[0]
r_eff = pearsonr(obs_eff, ppi_eff)[0]

fig_e, axes_e = plt.subplots(1, 2, figsize=(12.2, 5.4), constrained_layout=True)
lo, hi = float(min(obs_mean.min(), gen_mean.min())), float(max(obs_mean.max(), gen_mean.max()))
axes_e[0].scatter(obs_mean, gen_mean, alpha=0.4, s=12, c="steelblue", edgecolor="none")
axes_e[0].plot([lo, hi], [lo, hi], "r--", lw=1.0, alpha=0.65)
axes_e[0].set(xlabel="Real mean expression", ylabel="Predicted mean expression", title=f"Mean expression (Pearson r = {r_mean:.3f})")

lim_e = float(max(np.abs(obs_eff).max(), np.abs(ppi_eff).max())) * 1.05
axes_e[1].scatter(obs_eff, ppi_eff, alpha=0.4, s=12, c="coral", edgecolor="none")
axes_e[1].plot([-lim_e, lim_e], [-lim_e, lim_e], "r--", lw=1.0, alpha=0.65)
axes_e[1].axhline(0, color="grey", ls=":", lw=0.7)
axes_e[1].axvline(0, color="grey", ls=":", lw=0.7)
axes_e[1].set_xlim(-lim_e, lim_e)
axes_e[1].set_ylim(-lim_e, lim_e)
axes_e[1].set(xlabel="Real effect size (target - control)", ylabel="Predicted effect size", title=f"Effect size (Pearson r = {r_eff:.3f})")
fig_e.savefig("results/eval_mean_vs_effect_scatter.png", dpi=300, bbox_inches="tight")
plt.close(fig_e)
