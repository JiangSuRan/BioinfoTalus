## GSEA 多通路山峰图（Running ES 曲线叠加 + barcode 条带）
library(dplyr)
library(ggplot2)
library(patchwork)
library(enrichplot)

gsea_res <- readRDS("Data/gsea_result.rds")

selected_pathways_3B <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

pathway_colors <- c(
  "#E64B35",
  "#4DBBD5",
  "#00A087",
  "#3C5488",
  "#F39B7F",
  "#8491B4"
)

gene_names  <- names(gsea_res@geneList)
gsea_curve  <- data.frame()
barcode_dat <- data.frame()

for (pw in selected_pathways_3B) {
  tmp          <- enrichplot:::gsInfo(gsea_res, geneSetID = pw)
  tmp$id       <- pw
  gsea_curve   <- rbind(gsea_curve, tmp)
  barcode_dat  <- rbind(barcode_dat,
    tmp %>%
      dplyr::mutate(gene_name = gene_names) %>%
      dplyr::filter(position == 1)
  )
}

line_col <- setNames(pathway_colors, selected_pathways_3B)

label_3b <- gsea_curve %>%
  dplyr::group_by(id) %>%
  dplyr::slice(which.max(abs(runningScore))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = gsub("_", " ", gsub("^HALLMARK_", "", id))
  )

p3b_curve <- ggplot(gsea_curve, aes(x = x, y = runningScore)) +
  geom_line(aes(color = id), linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#555555", linewidth = 0.5) +
  geom_text(
    data        = label_3b,
    aes(x = x, y = runningScore, label = label, color = id),
    size        = 3,
    fontface    = "bold",
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = line_col,
    breaks = selected_pathways_3B,
    labels = setNames(
      gsub("_", " ", gsub("^HALLMARK_", "", selected_pathways_3B)),
      selected_pathways_3B
    ),
    name   = NULL
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(gsea_curve$x))) +
  labs(x = NULL, y = "Running Enrichment Score") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.line.x       = element_blank(),
    panel.grid        = element_blank(),
    legend.position   = c(0.85, 0.75),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    legend.background = element_rect(fill = "transparent"),
    plot.margin       = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm"),
    panel.border      = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

p3b_barcode <- ggplot(gsea_curve, aes(x = x, y = runningScore, color = id)) +
  geom_segment(
    data        = barcode_dat,
    aes(x = x, xend = x, y = 0, yend = 1),
    show.legend = FALSE,
    linewidth   = 0.3
  ) +
  scale_color_manual(values = line_col) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ id, ncol = 1) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_blank(),
    axis.line.x      = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.spacing.y  = unit(1, "mm"),
    plot.margin      = margin(t = 0, r = 0.2, b = 0, l = 0.2, unit = "cm"),
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

p3b_final <- patchwork::wrap_plots(p3b_curve, p3b_barcode,
                                   ncol = 1, heights = c(4, 2.6))

ggsave(filename = "gsea_mountain_plot.png", plot = p3b_final, width = 8, height = 6, units = "in", dpi = 300)
ggsave(filename = "gsea_mountain_plot.pdf", plot = p3b_final, width = 8, height = 6, units = "in")
