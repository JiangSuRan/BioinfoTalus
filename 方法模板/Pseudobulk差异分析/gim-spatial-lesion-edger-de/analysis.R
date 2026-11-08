## GIM-spatial 空间病理分区 Pseudobulk 聚合与 limma-voom 差异表达
library(edgeR)
library(limma)
library(tidyverse)

counts <- read_tsv("Data/pseudobulk_counts.tsv")
samples <- read_csv("Data/samples.csv")

counts_df <- counts %>%
  select(-(Gene)) %>%
  as.data.frame()
rownames(counts_df) <- counts %>% pull(Gene)

# 伪样本名按 array_LIMSdissection_group 拆分并关联患者表
pseudobulk_metadata <- tibble(pseudobulk = colnames(counts_df)) %>%
  separate(pseudobulk, c("array", "LIMS_dissection", "group"), sep = "_", remove = FALSE) %>%
  mutate(LIMS_dissection = as.double(LIMS_dissection)) %>%
  left_join(samples, by = "LIMS_dissection")

group <- factor(pseudobulk_metadata$group, levels = c("Base", "Pit", "Metaplasia"))
y <- DGEList(counts = counts_df, group = group, samples = pseudobulk_metadata %>%
  mutate(dissection = factor(LIMS_dissection)) %>% as.data.frame())
y$samples$array <- factor(pseudobulk_metadata$array.x)
y$samples$dissection <- factor(pseudobulk_metadata$LIMS_dissection)
y$samples$region <- factor(pseudobulk_metadata$region)

# 低表达过滤 + TMM 归一化
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

# 模型：分区组间比较 + 患者（dissection）与阵列批次协变量
dissection <- y$samples$dissection
region <- y$samples$region
design <- model.matrix(~0 + group + dissection + region)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
  MetaplasiavsBase = Metaplasia - Base,
  MetaplasiavsPit = Metaplasia - Pit,
  BasevsPit = Base - Pit,
  levels = colnames(design)
)

# voom 异方差校正后线性拟合与经验贝叶斯检验（treat 报告 lfc 过滤结果）
v <- voom(y, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
tfit <- treat(vfit, lfc = 0)
dt <- decideTests(tfit)

res <- topTable(efit, coef = "MetaplasiavsBase", number = Inf)
res_treat <- topTable(tfit, coef = "MetaplasiavsBase", number = Inf)

write.csv(res, "results/spatial_region_limma_results.csv")
write.csv(res_treat, "results/spatial_region_treat_results.csv")
write.csv(as.data.frame(dt), "results/spatial_region_decideTests.csv")
saveRDS(y, "results/edger_dgelist.rds")
