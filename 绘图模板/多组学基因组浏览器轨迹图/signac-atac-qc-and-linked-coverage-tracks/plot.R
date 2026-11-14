## Signac scATAC 质控联合散点/核小体周期条纹图与多组学连锁 Coverage 轨迹图
library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

atac_obj <- readRDS("Data/multiome_atac_linked.rds")
DefaultAssay(atac_obj) <- "peaks"

# 核小体信号（147 bp 单核小体与无核小体片段比值）用于区分开放染色质带纹质量
atac_obj$nucleosome_group <- ifelse(atac_obj$nucleosome_signal > 4, "NS > 4", "NS <= 4")

p_density <- DensityScatter(
  atac_obj,
  x = "nCount_peaks",
  y = "TSS.enrichment",
  log_x = TRUE,
  quantiles = TRUE
) +
  ggtitle("Peak Counts vs TSS Enrichment")

p_frag <- FragmentHistogram(
  object = atac_obj,
  group.by = "nucleosome_group",
  region = "chr1-1-2000000"
) +
  ggtitle("Fragment Length Periodicity")

p_qc <- p_density + p_frag
ggsave(filename = "results/signac_atac_qc_panel.pdf", plot = p_qc, width = 11, height = 4.8, units = "in")

# 联合展示伪批量开放信号轨道、基因表达小提琴、Call Peak 及 Peak-to-Gene 连锁弧线
highlight_gr <- StringToGRanges("chr1-1010000-1012000", sep = c("-", "-"))
p_cov <- CoveragePlot(
  object = atac_obj,
  region = "ISG15",
  features = "ISG15",
  expression.assay = "SCT",
  peaks = TRUE,
  links = TRUE,
  region.highlight = highlight_gr,
  extend.upstream = 5000,
  extend.downstream = 5000
)

ggsave(filename = "results/signac_linked_coverage_track.pdf", plot = p_cov, width = 10, height = 7, units = "in")
