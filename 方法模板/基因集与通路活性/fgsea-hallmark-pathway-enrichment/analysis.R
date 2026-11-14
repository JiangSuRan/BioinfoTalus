## GSVA / Welch t 检验治疗响应分组通路富集评分与统计比较
library(GSVA)
library(BiocParallel)
library(dplyr)
library(tibble)
library(jsonlite)

log_abundance <- read.table(
  "Data/bulk_logtpm.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
metadata <- read.table("Data/bulk_metadata.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
json_data <- jsonlite::fromJSON("Data/c2.cp.reactome.v2023.1.Hs.json")

IFNG_signature <- json_data[["REACTOME_INTERFERON_GAMMA_SIGNALING"]][["geneSymbols"]]
geneset_list <- list(IFNG_signature = IFNG_signature)

# GSVA 2.x 参数对象接口；连续 log-TPM 矩阵选用 Gaussian 核密度估计
score_param <- gsvaParam(
  exprData = as.matrix(log_abundance),
  geneSets = geneset_list,
  minSize = 1,
  maxSize = Inf,
  kcdf = "Gaussian",
  maxDiff = TRUE
)
score <- gsva(score_param, verbose = FALSE, BPPARAM = SerialParam(progressbar = FALSE))

pdat <- score %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Ident") %>%
  left_join(metadata, by = "Ident") %>%
  mutate(int_lab = paste(SampleTimePoint, Efficacy, sep = "_"))

# 治疗时点 × 疗效分组间通路得分差异（Welch t 检验 + BH 校正）
comparisons <- list(c("BL_", "nRCT_non-CR"), c("BL_", "nRCT_CR"), c("nRCT_non-CR", "nRCT_CR"))
stats_tab <- do.call(rbind, lapply(comparisons, function(cmp) {
  x <- pdat$IFNG_signature[pdat$int_lab == cmp[1]]
  y <- pdat$IFNG_signature[pdat$int_lab == cmp[2]]
  tt <- t.test(x, y)
  data.frame(
    comparison = paste(cmp[1], cmp[2], sep = "_vs_"),
    mean_group1 = mean(x),
    mean_group2 = mean(y),
    t_stat = unname(tt$statistic),
    p.value = tt$p.value
  )
}))
stats_tab$padj <- p.adjust(stats_tab$p.value, method = "BH")

write.csv(stats_tab, "results/signature_score_stats.csv", row.names = FALSE)
