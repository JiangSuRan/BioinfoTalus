## SMR / OPERA 多组学（GWAS + cis-eQTL + cis-caQTL）孟德尔随机化基因座共定位图与效应量散点图
source("Data/plot_OmicsSMR_xQTL_with_legend.r")
source("Data/plot_SMR.r")

pdf("results/smr_locus_and_effect_plot.pdf", width = 10, height = 8)

# 多层基因座关联轨迹：GWAS 与 eQTL/caQTL 逐层比对，突出共定位探针
SMRData_opera <- ReadomicSMRData("Data/smr_locus_opera.txt")
omicSMRLocusPlot(
  data = SMRData_opera,
  esmr_thresh = 0.00028,
  msmr_thresh = 0.001,
  esmr_heidi = 0.1,
  msmr_heidi = 0.1,
  trait_name = "As",
  window = 500,
  highlight = "chr12_56050848"
)

# 效应量散点回归：xQTL beta 对 GWAS beta 的方向一致性与解释力度
SMRData_smr <- ReadSMRData("Data/smr_locus_smr.txt")
SMREffectPlot(data = SMRData_smr, trait_name = "As")

dev.off()
