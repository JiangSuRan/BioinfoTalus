## tcrdist CDR3 距离矩阵构建阈值邻接网络与 Louvain 抗原特异性 TCR 聚落识别
import os

import networkx as nx
import pandas as pd
import community.community_louvain as community_louvain
from tcrdist.public import _neighbors_fixed_radius
from tcrdist.repertoire import TCRrep

data_path = "Data/"
clones_file = pd.read_csv(os.path.join(data_path, "tcell_tcr_meta_filter.csv"))
clones_file = clones_file.drop(columns=["Unnamed: 0"])

# 目标细胞类型子集；V/J 基因补 *01 等位后缀
targetcell = ["CD8 effector"]
df = clones_file[clones_file["ct_level_4"].isin(targetcell)].copy()
df["v_b_gene"] = df["v_b_gene"].apply(lambda x: f"{x}*01")
df["j_b_gene"] = df["j_b_gene"].apply(lambda x: f"{x}*01")

clone_cols = ["Condition", "cdr3_b_aa", "v_b_gene", "j_b_gene", "ct_level_4", "barcode", "Donor"]
edge_threshold = 24

tr = TCRrep(cell_df=df[clone_cols], chains=["beta"], organism="human")

# 固定半径近邻展开为带权边列表（权重 = 1 - dist/threshold）
network = []
for i, n in enumerate(_neighbors_fixed_radius(tr.pw_beta, edge_threshold)):
    for j in n:
        if i != j:
            network.append((
                i, j,
                tr.pw_beta[i, j],
                tr.clone_df["cdr3_b_aa"].iloc[i],
                tr.clone_df["cdr3_b_aa"].iloc[j],
                tr.clone_df["count"].iloc[i],
                tr.clone_df["count"].iloc[j],
                tr.clone_df["Condition"].iloc[i],
                tr.clone_df["Condition"].iloc[j],
                len(n) - 1,
            ))

cols = ["node_1", "node_2", "dist", "cdr3_b_aa_1", "cdr3_b_aa_2",
        "count_1", "count_2", "Condition_1", "Condition_2", "K_neighbors"]
df_net = pd.DataFrame(network, columns=cols)
df_net["weight"] = (edge_threshold - df_net["dist"]) / edge_threshold

# Louvain 社区识别并按规模重排簇编号
G = nx.from_pandas_edgelist(pd.DataFrame({"source": df_net["node_1"], "target": df_net["node_2"], "weight": df_net["weight"]}))
partition = community_louvain.best_partition(G, resolution=1)
partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
partition_reorder = {id: rank for id, rank in zip(partitions_by_cluster_size, range(len(partitions_by_cluster_size)))}
partition = {k: partition_reorder.get(v) for k, v in partition.items()}

df_net["community"] = df_net["node_1"].map(partition)
df_net.to_csv("results/tcrdist_network_edges.csv", index=False)

nodes = sorted(set(df_net["node_1"]) | set(df_net["node_2"]))
node_table = tr.clone_df.iloc[nodes][["cdr3_b_aa", "v_b_gene", "count", "Condition"]].copy()
node_table["community"] = pd.Series(partition)
node_table.to_csv("results/tcrdist_network_nodes.csv")
