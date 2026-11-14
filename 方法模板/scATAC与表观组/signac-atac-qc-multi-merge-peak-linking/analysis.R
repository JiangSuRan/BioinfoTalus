## Signac scATAC 多样本峰值合并定量、核小体/TSS 质控与 Peak-to-Gene 调控连锁
library(Seurat)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(dplyr)

peaks_a <- read.table("Data/sample_a_peaks.bed", col.names = c("chr", "start", "end"))
peaks_b <- read.table("Data/sample_b_peaks.bed", col.names = c("chr", "start", "end"))
meta_a <- read.csv("Data/sample_a_singlecell.csv", row.names = 1)
meta_b <- read.csv("Data/sample_b_singlecell.csv", row.names = 1)

# 将多个独立样本 Call Peak 区间取并集并合并重叠区段，构建跨样本可比的一致性 Peak 特征矩阵
gr_a <- makeGRangesFromDataFrame(peaks_a)
gr_b <- makeGRangesFromDataFrame(peaks_b)
combined_peaks <- reduce(c(gr_a, gr_b))
peak_widths <- width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

frag_a <- CreateFragmentObject("Data/sample_a_fragments.tsv.gz", cells = rownames(meta_a))
frag_b <- CreateFragmentObject("Data/sample_b_fragments.tsv.gz", cells = rownames(meta_b))

mat_a <- FeatureMatrix(fragments = frag_a, features = combined_peaks, cells = rownames(meta_a))
mat_b <- FeatureMatrix(fragments = frag_b, features = combined_peaks, cells = rownames(meta_b))

# 直接通过 ensembldb 提取标准染色体外显子 GRanges，避免 Signac::GetGRangesFromEnsDb 对 biovizBase 的额外依赖
annotations <- ensembldb::exons(
  EnsDb.Hsapiens.v86,
  columns = c("tx_id", "gene_name", "gene_id", "gene_biotype")
)
annotations <- annotations[annotations$gene_biotype %in% c("protein_coding", "lincRNA", "rRNA", "processed_transcript")]
annotations$type <- "exon"
annotations <- GenomeInfoDb::keepStandardChromosomes(annotations, pruning.mode = "coarse")
GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"

assay_a <- CreateChromatinAssay(counts = mat_a, fragments = frag_a, annotation = annotations)
assay_b <- CreateChromatinAssay(counts = mat_b, fragments = frag_b, annotation = annotations)

obj_a <- CreateSeuratObject(counts = assay_a, assay = "peaks", meta.data = meta_a)
obj_b <- CreateSeuratObject(counts = assay_b, assay = "peaks", meta.data = meta_b)
obj_a$sample_id <- "SampleA"
obj_b$sample_id <- "SampleB"

atac <- merge(obj_a, y = obj_b, add.cell.ids = c("SampleA", "SampleB"))
atac <- NucleosomeSignal(atac)
atac <- TSSEnrichment(atac, fast = FALSE)
atac$blacklist_fraction <- FractionCountsInRegion(atac, assay = "peaks", regions = blacklist_hg38_unified)

atac <- subset(
  atac,
  subset = nCount_peaks > 3000 & nCount_peaks < 30000 & nucleosome_signal < 4 & TSS.enrichment > 2 & blacklist_fraction < 0.05
)

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30)
atac <- FindNeighbors(atac, reduction = "lsi", dims = 2:30)
atac <- FindClusters(atac, algorithm = 3, resolution = 0.5)

# 基于启动子及基因体开放片段估算 GeneActivity，便于与配对表达层共同进行调控联动推断
gene_activities <- GeneActivity(atac)
atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activities)
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)

saveRDS(atac, "results/signac_qc_merged_atac.rds")
