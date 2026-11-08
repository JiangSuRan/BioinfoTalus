## CellOracle 基因调控网络（GRN）拟合与转录因子虚拟敲除（In-silico KO）分化轨迹模拟
import celloracle as co

oracle = co.load_hdf5("Data/oracle_object.celloracle")
links = co.load_hdf5("Data/links_object.celloracle")

# 过滤边并按簇特异 TF 字典重拟合岭回归模拟模型
links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

# 虚拟敲除：目标 TF 表达置 0，经 3 步 GRN 信号传播得到扰动后全基因表达
goi = "Gata1"
oracle.simulate_shift(perturb_condition={goi: 0.0}, n_propagation=3)

# 随机图近似的马尔可夫转移概率与嵌入空间偏移向量
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

oracle.to_hdf5("results/oracle_gata1_ko.celloracle")
