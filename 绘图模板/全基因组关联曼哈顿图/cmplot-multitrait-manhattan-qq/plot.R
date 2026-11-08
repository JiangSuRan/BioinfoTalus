## CMplot 多表型全基因组关联分析（GWAS）背景降采样矩形曼哈顿图与多轨迹分位数 QQ 图
library(CMplot)
library(data.table)
library(dplyr)

traits <- c("HDL", "LDL", "TC", "TG")

# 逐表型读取 GEMMA 结果并合并 p 值列
gwas_df <- fread(paste0("Data/gwas/", traits[1], "_378.mlma"), sep = "\t")
gwas_df <- gwas_df[, c("SNP", "Chr", "bp")]
colnames(gwas_df) <- c("SNP", "Chromosome", "Position")
for (tr in traits) {
  temp_gwas <- fread(paste0("Data/gwas/", tr, "_378.mlma"), sep = "\t")
  gwas_df[, (tr) := temp_gwas$p]
}
gwas_df <- as.data.frame(gwas_df)

# 全基因组显著位点全保留，非显著背景降采样避免千万点拖垮矢量输出
sig_df <- gwas_df %>% filter(if_any(all_of(traits), ~ .x < 5e-8))
no_sig_df <- gwas_df %>% filter(if_all(all_of(traits), ~ .x >= 5e-8))
df_plot <- rbind(sig_df, no_sig_df %>% sample_frac(0.02))

CMplot(
  df_plot,
  plot.type = "m",
  col = c("gray50"),
  multracks = TRUE,
  threshold = c(5e-8, 5e-8 / 36),
  threshold.lty = 1,
  threshold.lwd = c(1, 1),
  threshold.col = c("black", "grey"),
  amplify = TRUE,
  chr.den.col = NULL,
  signal.col = c("#011627", "#F71735", "#41EAD4", "#FF9F1C"),
  signal.cex = 1,
  file = "pdf",
  file.name = "results/cima_multitrait_manhattan",
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  points.alpha = 225,
  legend.ncol = 3,
  legend.pos = "middle"
)

CMplot(
  df_plot,
  plot.type = "q",
  col = c("#011627", "#F71735", "#41EAD4", "#FF9F1C"),
  multracks = TRUE,
  ylab.pos = 2,
  signal.pch = c(19, 6, 4),
  signal.cex = 1.2,
  signal.col = "red",
  conf.int = TRUE,
  box = FALSE,
  axis.cex = 1,
  file = "pdf",
  file.name = "results/cima_multitrait_qq",
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 8,
  height = 8
)

write.csv(sig_df, "results/cima_gwas_significant_snps.csv", row.names = FALSE)
