## SAHMI 单细胞宿主-微生物组双端条形码降噪过滤与寄主细胞群整合分析
library(Seurat)
library(data.table)
library(dplyr)

sc_obj <- readRDS("Data/host_seurat.rds")
species_counts <- fread("Data/sahmi_cell_species_counts.tsv")
kmer_metrics <- fread("Data/sahmi_kmer_consistency.tsv")

# 过滤假阳性环境污染微生物：要求跨细胞重复检出、唯一 k-mer 数与 read 数保持单调正相关
valid_taxa <- kmer_metrics %>%
  filter(
    total_reads >= 50,
    total_unique_kmers >= 100,
    n_cells_detected >= 5,
    spearman_kmer_read >= 0.6,
    pvalue_kmer_read < 0.05
  ) %>%
  pull(taxon_name)

filtered_counts <- species_counts %>%
  filter(taxon_name %in% valid_taxa, cell_barcode %in% colnames(sc_obj))

cell_microbe_summary <- filtered_counts %>%
  group_by(cell_barcode) %>%
  summarise(
    microbe_umi = sum(umi_count),
    dominant_taxon = taxon_name[which.max(umi_count)],
    .groups = "drop"
  )

meta_add <- data.frame(
  row.names = colnames(sc_obj),
  microbe_umi = rep(0, ncol(sc_obj)),
  dominant_taxon = rep("Uncolonized", ncol(sc_obj)),
  stringsAsFactors = FALSE
)
matched <- match(cell_microbe_summary$cell_barcode, rownames(meta_add))
meta_add$microbe_umi[matched] <- cell_microbe_summary$microbe_umi
meta_add$dominant_taxon[matched] <- cell_microbe_summary$dominant_taxon
meta_add$microbe_status <- ifelse(meta_add$microbe_umi > 0, "Positive", "Negative")

sc_obj <- AddMetaData(sc_obj, metadata = meta_add)
Idents(sc_obj) <- "microbe_status"
host_response_deg <- FindMarkers(sc_obj, ident.1 = "Positive", ident.2 = "Negative", min.pct = 0.1)

saveRDS(sc_obj, "results/sahmi_host_microbe_seurat.rds")
write.csv(host_response_deg, "results/sahmi_colonized_vs_uncolonized_deg.csv", row.names = TRUE)
