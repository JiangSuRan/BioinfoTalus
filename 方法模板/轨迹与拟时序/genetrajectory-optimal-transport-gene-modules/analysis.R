## GeneTrajectory 基于细胞图最优传输（Wasserstein）距离推断基因时序动态模块
library(Seurat)
library(GeneTrajectory)

sc_obj <- readRDS("Data/seurat_preprocessed.rds")
selected_genes <- readLines("Data/variable_trajectory_genes.txt")

# 构建扩散映射（Diffusion Map）捕捉连续细胞流形结构，并计算细胞 kNN 最短图距离
sc_obj <- RunDM(sc_obj)
cell_graph_dist <- GetGraphDistance(sc_obj, K = 10)

# 将大规模细胞图粗粒化（CoarseGrain）为代表性元节点，以降低基因间最优传输距离计算复杂度
cg_out <- CoarseGrain(
  sc_obj,
  cell_graph_dist,
  genes = intersect(selected_genes, rownames(sc_obj)),
  N = 500
)
gene_ot_dist <- CalGeneDistByOT(cg_out)

# 构建基因-基因亲和图谱与扩散嵌入，提取有序共调控基因轨迹（Bin）并映射回单细胞打分
gene_embed <- GetGeneEmbedding(gene_ot_dist, K = 5)
gene_traj <- ExtractGeneTrajectory(gene_embed, gene_ot_dist, N = 3, t.list = c(5, 5, 5))
sc_obj <- AddGeneBinScore(sc_obj, gene_traj, N.bins = 5, trajectories = 1:3)

saveRDS(sc_obj, "results/genetrajectory_scored_seurat.rds")
saveRDS(gene_traj, "results/genetrajectory_modules.rds")
