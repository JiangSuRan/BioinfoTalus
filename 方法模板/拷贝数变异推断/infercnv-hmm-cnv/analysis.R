## inferCNV 单细胞/空间转录组拷贝数变异推断（HMM 降噪模式）
library(infercnv)

counts <- readRDS("Data/raw_counts_matrix.rds")

# ref_group_names 指定正常参考群，以其表达基线校准目标群（肿瘤/未知）的拷贝数推断
infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts,
  gene_order_file = "Data/gene_order.tsv",
  annotations_file = "Data/cell_annotations.tsv",
  delim = "\t",
  ref_group_names = c("Base", "Pit")
)

# HMM 将校正后的表达离散化为 CNV 状态，denoise + noise_filter 压制小幅波动假阳性
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "results/infercnv_out",
  cluster_by_groups = FALSE,
  scale_data = TRUE,
  num_threads = 8,
  denoise = TRUE,
  noise_filter = 0.2,
  HMM = TRUE
)
