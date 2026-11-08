## CrossTalkeR 跨条件配体-受体通讯网络比较与差异中心性分析
library(Seurat)
library(liana)
library(CrossTalkeR)
library(dplyr)

adata_seurat <- readRDS("Data/athero_seurat.rds")
adata_seurat <- subset(adata_seurat, subset = (Main_Cell_Types != "Osteoblastic cells"))

data <- as.SingleCellExperiment(adata_seurat)

# 按条件拆分 counts 与 metadata，分别构建 Seurat 对象
gp <- assays(data)[["counts"]][, colData(data)$Disease_State == "Control"]
pp <- assays(data)[["counts"]][, colData(data)$Disease_State == "Atherosclerosis"]

metadata_gp <- tibble::tibble(
  rownames(colData(data))[colData(data)$Disease_State == "Control"],
  as.character(colData(data)$Main_Cell_Types[colData(data)$Disease_State == "Control"])
)
metadata_pp <- tibble::tibble(
  rownames(colData(data))[colData(data)$Disease_State == "Atherosclerosis"],
  as.character(colData(data)$Main_Cell_Types[colData(data)$Disease_State == "Atherosclerosis"])
)

good <- CreateSeuratObject(counts = gp)
Idents(good) <- metadata_gp$`as.character(...)`

poor <- CreateSeuratObject(counts = pp)
Idents(poor) <- metadata_pp$`as.character(...)`

# 两条件分别运行 LIANA 并取 cellphonedb 结果为 CrossTalkeR 输入格式
good <- NormalizeData(good, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test <- liana_wrap(good)

poor <- NormalizeData(poor, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test2 <- liana_wrap(poor)

control <- liana_test$cellphonedb %>%
  filter(pvalue <= 0.05) %>%
  mutate(gene_A = gsub("_", "*", ligand.complex)) %>%
  mutate(gene_B = gsub("_", "*", receptor.complex)) %>%
  mutate(source = gsub("_", "", source)) %>%
  mutate(target = gsub("_", "", target)) %>%
  mutate(type_gene_A = rep("Ligand", length(.data$ligand))) %>%
  mutate(type_gene_B = rep("Receptor", length(.data$receptor)))

acm <- liana_test2$cellphonedb %>%
  filter(pvalue <= 0.05) %>%
  mutate(gene_A = gsub("_", "*", ligand.complex)) %>%
  mutate(gene_B = gsub("_", "*", receptor.complex)) %>%
  mutate(source = gsub("_", "", source)) %>%
  mutate(target = gsub("_", "", target)) %>%
  mutate(type_gene_A = rep("Ligand", length(.data$ligand))) %>%
  mutate(type_gene_B = rep("Receptor", length(.data$receptor)))

l <- list()
l[["Control"]] <- control
l[["Athero"]] <- acm

# 生成跨条件比较 HTML 报告（PageRank/介数中心性差异与差异通讯边）
dir.create("results/crosstalker", showWarnings = FALSE)
data_dis <- generate_report(
  lrpaths = l,
  genes = NULL,
  out_path = "results/crosstalker/",
  threshold = 0,
  out_file = "result.html",
  output_fmt = "html_document",
  report = TRUE,
  sel_columns = c("source", "target", "gene_A", "gene_B", "type_gene_A", "type_gene_B", "lr.mean")
)
