## TCGA STAR counts 整理
library(tidyverse)
library(jsonlite)

expr_dir <- "Data/TCGA/RNA seq/counts files"
meta_json <- "Data/TCGA/RNA seq/metadata.cart.2025-09-22.json"
out_dir <- "Data/TCGA_expression_output"
keep_sample_types <- c("01" = "Primary Tumor", "11" = "Solid Tissue Normal")
derive_symbol_matrix <- TRUE
write_text <- FALSE

stopifnot(dir.exists(expr_dir))
files <- list.files(expr_dir, pattern = "gene_counts\\.tsv(\\.gz)?$",
                    recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) stop("未找到 STAR gene_counts.tsv 或 .tsv.gz 文件")

# 先读原始行定位表头，避免注释行与数据列数不同影响行号
read_star_counts <- function(path) {
  lines <- readLines(path, n = 10, warn = FALSE)
  header_line <- grep("^gene_id\t", lines)
  if (length(header_line) != 1) stop("前 10 行未找到唯一 gene_id 表头：", path)
  dat <- read.delim(path, skip = header_line - 1, check.names = FALSE)
  required <- c("gene_id", "gene_name", "gene_type", "unstranded", "tpm_unstranded")
  if (!all(required %in% names(dat))) stop("STAR 文件缺少必需字段：", path)
  dat <- dat[!dat$gene_id %in% c("N_unmapped", "N_multimapping",
                                "N_noFeature", "N_ambiguous"), ]
  if (anyNA(dat$gene_id) || any(dat$gene_id == "") || anyDuplicated(dat$gene_id)) {
    stop("gene_id 缺失或重复：", path)
  }
  dat
}

ref <- read_star_counts(files[1])
gene_annotation <- ref %>% dplyr::select(gene_id, gene_name, gene_type)
gene_ids <- ref$gene_id

# 矩阵每列对应一个文件，样本表每行对应一个文件
per_file <- tibble(file_dir = basename(dirname(files)),
                   file_name = basename(files), lib_size = NA_real_)
if (anyDuplicated(per_file$file_dir)) stop("每个文件需位于独立目录，目录名用于矩阵列键")
counts_all <- matrix(NA_real_, nrow = length(gene_ids), ncol = length(files),
                     dimnames = list(gene_ids, per_file$file_dir))
tpm_all <- counts_all
for (i in seq_along(files)) {
  dat <- read_star_counts(files[i])
  if (!setequal(dat$gene_id, gene_ids)) stop("基因集不一致：", files[i])
  dat <- dat[match(gene_ids, dat$gene_id), ]
  if (!identical(dat$gene_name, ref$gene_name) ||
      !identical(dat$gene_type, ref$gene_type)) stop("基因注释不一致：", files[i])
  counts <- as.numeric(dat$unstranded)
  tpm <- as.numeric(dat$tpm_unstranded)
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stop("counts 不是非负整数：", files[i])
  }
  if (any(!is.finite(tpm)) || any(tpm < 0)) stop("TPM 缺失或不是非负有限值：", files[i])
  counts_all[, i] <- counts
  tpm_all[, i] <- tpm
  per_file$lib_size[i] <- sum(counts)
}

# metadata 和 API 返回使用同一套字段展开
map_entities <- function(hits) {
  rows <- list()
  for (h in hits) {
    for (e in h$associated_entities) {
      rows[[length(rows) + 1]] <- tibble(file_id = h$file_id,
                                        file_name = h$file_name,
                                        entity_barcode = e$entity_submitter_id)
    }
  }
  bind_rows(rows)
}

if (file.exists(meta_json)) {
  hits <- fromJSON(meta_json, simplifyVector = FALSE)
  file_map <- map_entities(hits)
} else {
  if (!requireNamespace("httr", quietly = TRUE)) stop("缺少 metadata；API 查询需安装 httr")
  uuid_chunks <- split(per_file$file_dir, ceiling(seq_len(nrow(per_file)) / 500))
  map_parts <- vector("list", length(uuid_chunks))
  for (i in seq_along(uuid_chunks)) {
    uuids <- uuid_chunks[[i]]
    filters <- toJSON(list(op = "in", content = list(field = "files.file_id",
                                                     value = as.list(uuids))), auto_unbox = TRUE)
    resp <- httr::POST("https://api.gdc.cancer.gov/files",
                       body = list(filters = filters,
                                   fields = "file_id,file_name,associated_entities.entity_submitter_id",
                                   size = length(uuids)),
                       encode = "form", httr::timeout(120))
    httr::stop_for_status(resp)
    hits <- fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                     simplifyVector = FALSE)$data$hits
    map_parts[[i]] <- map_entities(hits)
  }
  file_map <- bind_rows(map_parts)
}
if (!all(c("file_id", "file_name", "entity_barcode") %in% names(file_map))) {
  stop("metadata 缺少文件或 entity 映射字段")
}

uuid_map <- file_map %>% dplyr::select(file_id, entity_barcode)
sample_meta <- left_join(per_file, uuid_map, by = c("file_dir" = "file_id"))
if (anyNA(sample_meta$entity_barcode) && all(per_file$file_name %in% file_map$file_name)) {
  name_keys <- file_map %>% distinct(file_name, file_id)
  if (anyDuplicated(name_keys$file_name)) stop("文件名对应多个 UUID，不能按文件名回退匹配")
  name_map <- file_map %>% dplyr::select(file_name, entity_barcode)
  sample_meta <- left_join(per_file, name_map, by = "file_name")
}
if (anyNA(sample_meta$entity_barcode)) stop("部分文件无法映射 barcode，请检查 metadata")
if (any(!grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}[A-Za-z]($|-)",
               sample_meta$entity_barcode))) stop("entity 标识不是预期 TCGA barcode")

# 保留完整 entity ID，另取 16 位样本条码以落实原模板的去重规则
sample_meta <- sample_meta %>%
  mutate(sample_barcode = substr(entity_barcode, 1, 16),
         patient_barcode = substr(entity_barcode, 1, 12),
         sample_type_code = substr(entity_barcode, 14, 15),
         vial = substr(entity_barcode, 16, 16))
sample_meta$sample_type <- unname(keep_sample_types[sample_meta$sample_type_code])
sample_meta$sample_type[is.na(sample_meta$sample_type)] <- "Other"

# 保留原模板决策：多 entity 取首条、同样本取文库最大、同类取最小 vial
sample_meta$multi_entity <- duplicated(sample_meta$file_dir)
sample_meta <- sample_meta %>% arrange(sample_barcode, desc(lib_size), file_name)
# 已排除的 entity 不参与后续择优，避免其 vial 或文库大小排挤有效样本
eligible <- which(!sample_meta$multi_entity &
                    sample_meta$sample_type_code %in% names(keep_sample_types))
sample_meta$dup_file <- FALSE
sample_meta$dup_file[eligible] <- duplicated(sample_meta$sample_barcode[eligible])
vial_choices <- sample_meta %>%
  filter(!multi_entity, !dup_file, sample_type_code %in% names(keep_sample_types)) %>%
  group_by(patient_barcode, sample_type_code) %>%
  summarise(vial_keep = min(vial), .groups = "drop")
sample_meta <- left_join(sample_meta, vial_choices, by = c("patient_barcode", "sample_type_code"))
sample_meta$keep <- !sample_meta$multi_entity & !sample_meta$dup_file &
  sample_meta$vial == sample_meta$vial_keep &
  sample_meta$sample_type_code %in% names(keep_sample_types)
sample_meta <- sample_meta %>%
  mutate(removal_reason = case_when(
    multi_entity ~ "同一文件多个 entity，保留第一个",
    dup_file ~ "同一 16 位样本多个文件，保留文库最大者",
    !sample_type_code %in% names(keep_sample_types) ~ "sample_type 不在保留列表",
    vial != vial_keep ~ paste0("同类样本多 vial，保留 ", vial_keep),
    TRUE ~ NA_character_))
kept <- sample_meta %>% filter(keep)
if (nrow(kept) == 0 || anyDuplicated(kept$sample_barcode)) stop("保留样本为空或条码重复")

idx <- match(kept$file_dir, colnames(counts_all))
counts_kept <- counts_all[, idx, drop = FALSE]
tpm_kept <- tpm_all[, idx, drop = FALSE]
colnames(counts_kept) <- kept$sample_barcode
colnames(tpm_kept) <- kept$sample_barcode
pc <- gene_annotation$gene_type == "protein_coding"
counts_protein_coding <- counts_kept[pc, , drop = FALSE]
tpm_protein_coding <- tpm_kept[pc, , drop = FALSE]

if (derive_symbol_matrix) {
  symbols <- gene_annotation$gene_name[pc]
  if (anyNA(symbols) || any(symbols == "")) stop("symbol 缺失，不能生成 symbol 矩阵")
  counts_pc_symbol <- rowsum(counts_protein_coding, group = symbols, reorder = FALSE)
}

#+ Check Point：矩阵按 gene_id 对齐，列为按原选择规则保留的 16 位样本条码
message("保留样本：", ncol(counts_kept), "；患者：", n_distinct(kept$patient_barcode))
print(table(kept$sample_type))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(counts_kept, file.path(out_dir, "counts_all.rds"))
saveRDS(tpm_kept, file.path(out_dir, "tpm_all.rds"))
saveRDS(counts_protein_coding, file.path(out_dir, "counts_protein_coding.rds"))
saveRDS(tpm_protein_coding, file.path(out_dir, "tpm_protein_coding.rds"))
write_tsv(gene_annotation, file.path(out_dir, "gene_annotation.tsv"))
write_tsv(sample_meta, file.path(out_dir, "sample_metadata.tsv"))
if (derive_symbol_matrix) saveRDS(counts_pc_symbol, file.path(out_dir, "counts_protein_coding_symbol.rds"))
if (write_text) write.csv(counts_protein_coding, file.path(out_dir, "counts_protein_coding.csv"))
qc_lines <- c(paste0("STAR 文件数: ", length(files)),
              paste0("基因数: ", length(gene_ids), "；protein_coding: ", sum(pc)),
              paste0("保留样本: ", nrow(kept), "；患者: ", n_distinct(kept$patient_barcode)),
              "样本删减原因:", capture.output(print(table(sample_meta$removal_reason))),
              "保留样本类型:", capture.output(print(table(kept$sample_type))))
writeLines(qc_lines, file.path(out_dir, "qc_report.txt"))
