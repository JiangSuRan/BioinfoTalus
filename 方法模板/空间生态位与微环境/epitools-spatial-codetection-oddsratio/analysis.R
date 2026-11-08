## epitools Xenium SNV 基因型与肿瘤克隆共检出的列联表构建与比值比检验
library(dplyr)
library(epitools)
library(effectsize)

variant_key <- "TPD52L1-p-A21---ALT-T"
wt_key <- "TPD52L1-p-A21---WT"
assay <- "Xenium.with.snvs"

obj <- readRDS("Data/xenium_snv_seurat.rds")

# 逐细胞提取变体与野生型探针 UMI 并二值化（>0 计 ALT/WT）
counts_df <- t(as.matrix(GetAssayData(object = obj, assay = assay, layer = "counts")))
counts_df <- as.data.frame(counts_df)
obj@meta.data[[paste0(variant_key, "_", assay, "_count")]] <- counts_df[, variant_key]
obj@meta.data[[paste0(wt_key, "_", assay, "_count")]] <- counts_df[, wt_key]
call_col <- paste0(variant_key, "_", assay, "_call")
obj@meta.data[[call_col]] <- NA
obj@meta.data[[call_col]][obj@meta.data[[paste0(wt_key, "_", assay, "_count")]] > 0] <- "WT"
obj@meta.data[[call_col]][obj@meta.data[[paste0(variant_key, "_", assay, "_count")]] > 0] <- "ALT"

# 肿瘤克隆与正常标签
obj$tumor_cell_call <- obj$neoplasm_normal_unknown
levels(obj$tumor_cell_call) <- c("Normal", "Cancer_clone_1", "Cancer_clone_2", "low_quality/unknown")

# 克隆 × 基因型 2x2 列联表（按细胞求和探针 UMI）
data_1 <- data.frame(
  ALT = c(0, 0),
  WT = c(0, 0),
  row.names = c("Cancer_clone_2", "Normal")
)
data_1["Normal", "WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"), paste0(wt_key, "_", assay, "_count")])
data_1["Normal", "ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Normal"), paste0(variant_key, "_", assay, "_count")])
data_1["Cancer_clone_2", "WT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"), paste0(wt_key, "_", assay, "_count")])
data_1["Cancer_clone_2", "ALT"] <- sum(obj@meta.data[(obj@meta.data$tumor_cell_call == "Cancer_clone_2"), paste0(variant_key, "_", assay, "_count")])

# 比值比 + Fisher 置信区间，并做模拟 p 值卡方交叉验证
oddsratio.out <- oddsratio(as.matrix(data_1), method = "fisher")
variant_chisq_test <- chisq.test(data_1, simulate.p.value = TRUE, correct = FALSE)

res <- list(
  oddsratio = oddsratio.out,
  chisq_sim = variant_chisq_test
)
saveRDS(res, "results/codetection_oddsratio_res.rds")
write.csv(as.data.frame(oddsratio.out$measure), "results/codetection_oddsratio.csv")
