## GSEApy 富集气泡图、发散 NES 条形图与通路 Jaccard 重叠网络图
import gseapy as gp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

enr_res = pd.read_csv("Data/gseapy_enrichment_results.csv")

# 气泡图展示 FDR 与基因比例，enrichment_map 提取通路间重叠基因网络以识别冗余功能簇
fig, ax = plt.subplots(figsize=(7, 8))
gp.dotplot(
    enr_res,
    column="Adjusted P-value",
    x="Gene_set",
    size=8,
    top_term=15,
    cmap="viridis_r",
    ax=ax,
    show_ring=True,
)
fig.savefig("results/gseapy_enrichment_dotplot.pdf", dpi=300, bbox_inches="tight")

nodes, edges = gp.enrichment_map(enr_res.head(25))
G = nx.from_pandas_edgelist(edges, source="src_idx", target="targ_idx", edge_attr=["jaccard_coef", "overlap_coef"])
pos = nx.kamada_kawai_layout(G)

fig_net, ax_net = plt.subplots(figsize=(9, 8))
nx.draw_networkx_edges(
    G,
    pos=pos,
    width=[d["jaccard_coef"] * 10 for _, _, d in G.edges(data=True)],
    edge_color="#B8C4D4",
    alpha=0.75,
    ax=ax_net,
)
nx.draw_networkx_nodes(
    G,
    pos=pos,
    nodelist=list(nodes.index),
    node_size=nodes["Hits_ratio"] * 900,
    node_color=-np.log10(nodes["Adjusted P-value"].clip(lower=1e-300)),
    cmap="YlOrRd",
    ax=ax_net,
)
nx.draw_networkx_labels(
    G,
    pos=pos,
    labels=nodes["Term"].to_dict(),
    font_size=7,
    ax=ax_net,
)
ax_net.set_axis_off()
fig_net.savefig("results/gseapy_enrichment_map_network.pdf", dpi=300, bbox_inches="tight")
