## 候选基因定向 GRN 双环中心局部网络图与双向注意力关联热图
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

neigh_df = pd.read_csv("Data/scprint_candidate_neighbors.csv")
mat_df = pd.read_csv("Data/scprint_grn_matrix.csv", index_col=0)

cand_names = list(dict.fromkeys(neigh_df["candidate"].astype(str)))
top_df = neigh_df.sort_values("weight", ascending=False).head(36)

G = nx.Graph()
for _, row in top_df.iterrows():
    G.add_edge(str(row["candidate"]), str(row["neighbor"]), weight=float(row["weight"]))

cand_set = set(cand_names)
cands = [c for c in cand_names if c in G]
nbrs = [n for n in G.nodes() if n not in cand_set]

# 候选基因排入内环，邻居按所连候选角度的加权圆均值排序后交替排入两层外环以防标签重叠
pos, cand_ang = {}, {}
if len(cands) == 1:
    pos[cands[0]] = np.array([0.0, 0.0])
    cand_ang[cands[0]] = np.pi / 2
else:
    for i, c in enumerate(cands):
        a = 2 * np.pi * i / len(cands) + np.pi / 2
        cand_ang[c] = a
        pos[c] = np.array([1.5 * np.cos(a), 1.5 * np.sin(a)])

pref, shared = {}, {}
for n in nbrs:
    vx = vy = 0.0
    deg = 0
    for c in cands:
        if G.has_edge(n, c):
            w = G[n][c]["weight"]
            vx += w * np.cos(cand_ang[c])
            vy += w * np.sin(cand_ang[c])
            deg += 1
    pref[n] = np.arctan2(vy, vx) if (vx or vy) else 0.0
    shared[n] = deg >= 2

nbrs_sorted = sorted(nbrs, key=lambda n: pref[n])
for j, n in enumerate(nbrs_sorted):
    a = 2 * np.pi * j / len(nbrs_sorted) + np.pi / 2
    r = 3.0 + (0.55 if (j % 2) else 0.0)
    pos[n] = np.array([r * np.cos(a), r * np.sin(a)])

node_list = list(G.nodes())
cand_sizes = {c: 1700 + 220 * max(0, len(c) - 3) for c in cands}
node_colors = ["#e74c3c" if nd in cand_set else "#5dade2" for nd in node_list]
node_sizes = [cand_sizes[nd] if nd in cand_set else 470 for nd in node_list]
max_w = max(G[u][v]["weight"] for u, v in G.edges())

fig, ax = plt.subplots(figsize=(11, 9.5), constrained_layout=True)
for u, v in G.edges():
    w = G[u][v]["weight"]
    ax.plot(
        [pos[u][0], pos[v][0]],
        [pos[u][1], pos[v][1]],
        color="#aab7c4",
        lw=0.7 + 4.2 * (w / max_w),
        alpha=0.22 + 0.6 * (w / max_w),
        zorder=1,
        solid_capstyle="round",
    )

nx.draw_networkx_nodes(
    G, pos, nodelist=node_list, node_color=node_colors,
    node_size=node_sizes, alpha=0.95, edgecolors="white", linewidths=1.3, ax=ax,
)
for c in cands:
    ax.text(
        pos[c][0], pos[c][1], c,
        fontsize=(12 if len(c) <= 4 else 10), fontweight="bold",
        color="white", ha="center", va="center", zorder=5,
    )
for n in nbrs:
    x, y = pos[n]
    rr = np.hypot(x, y)
    scale = 1.0 + 0.34 / max(rr, 1e-6)
    ax.text(
        x * scale, y * scale, n,
        fontsize=9, ha="center", va="center", color="#1b2631", zorder=4,
        fontweight=("bold" if shared[n] else "normal"),
    )

legend = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#e74c3c", markersize=13, label="Candidate gene (focal)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#5dade2", markersize=10, label="Associated neighbor gene"),
    Line2D([0], [0], color="#aab7c4", lw=3, label="Association (thicker = stronger)"),
]
ax.legend(handles=legend, loc="upper left", fontsize=10, framealpha=0.92)
ax.set_title(f"Candidate-centered GRN ({G.number_of_edges()} top associations)", fontsize=14)
ax.axis("off")
ax.set_aspect("equal")
ax.margins(0.14)
fig.savefig("results/candidate_local_network.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# 有符号 GRN 注意力关联矩阵热图（按非零权重第 95 百分位对称截断）
mat = mat_df.to_numpy(dtype=float)
genes = mat_df.index.astype(str).tolist()
nz = np.abs(mat[mat != 0])
vmax = float(np.percentile(nz, 95))

fig_hm, ax_hm = plt.subplots(figsize=(8.2, 7.0), constrained_layout=True)
im = ax_hm.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
ax_hm.set_xticks(range(len(genes)))
ax_hm.set_yticks(range(len(genes)))
ax_hm.set_xticklabels(genes, rotation=90, fontsize=8.5)
ax_hm.set_yticklabels(genes, fontsize=8.5)
ax_hm.set_xticks(np.arange(-0.5, len(genes), 1), minor=True)
ax_hm.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
ax_hm.grid(which="minor", color="lightgray", linewidth=0.5)
ax_hm.tick_params(which="minor", length=0)
ax_hm.set(xlabel="Target gene", ylabel="Regulator gene", title="Gene association matrix (attention weight)")
fig_hm.colorbar(im, ax=ax_hm, shrink=0.85, label="Association strength")
fig_hm.savefig("results/grn_association_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig_hm)
