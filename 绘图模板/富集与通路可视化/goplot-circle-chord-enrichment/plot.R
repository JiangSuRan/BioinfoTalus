## GO 通路富集 z-score 圆环图与基因-通路关联弦图
library(GOplot)
library(stringr)
library(dplyr)

enrich_obj <- readRDS("Data/enrich_go_bp.rds")
deg_df <- readRDS("Data/deg_table.rds")

go_df <- enrich_obj@result %>%
  slice_head(n = 3) %>%
  transmute(
    Category = "BP",
    ID = ID,
    Term = Description,
    Genes = str_replace_all(geneID, "/", ","),
    adj_pval = p.adjust
  )

term_genes <- strsplit(go_df$Genes, ",", fixed = TRUE) %>% unlist() %>% unique()
valid_genes <- intersect(term_genes, rownames(deg_df))
gene_df <- data.frame(ID = valid_genes, logFC = deg_df[valid_genes, "avg_log2FC"])

circ_dat <- circle_dat(go_df, gene_df)
chord_mat <- chord_dat(circ_dat, gene_df)

pdf("results/goplot_circle_chord.pdf", width = 10, height = 9)
print(GOCircle(circ_dat, nsub = nrow(go_df)))
print(GOChord(chord_mat, gene.order = "logFC"))
dev.off()
