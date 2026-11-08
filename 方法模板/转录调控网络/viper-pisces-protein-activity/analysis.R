## VIPER / PISCES 基于 ARACNe 调控网络的蛋白质/转录因子活性推断
library(viper)

# 整合 ARACNe bootstrap 网络：逐轮累积每条边的方向符号与互等信息
collapse.key <- "--aracne-edge--"
wd <- "Data/aracne_bootstrap"
fns <- list.files(path = wd, pattern = "bootstrapNetwork")
res <- resmi <- edges <- NULL
for (fn1 in fns) {
  tmp <- strsplit(readLines(file.path(wd, fn1)), "\t")
  mi <- as.numeric(sapply(tmp, function(x) x[3]))
  names(mi) <- sapply(tmp, function(x) paste(x[1:2], collapse = collapse.key))
  edges <- c(edges, length(mi))
  ed <- unique(c(names(res), names(mi)))
  pos1 <- match(ed, names(res))
  pos1na <- which(is.na(pos1))
  pos2 <- match(ed, names(mi))
  pos2na <- which(is.na(pos2))
  res <- res[pos1]
  resmi <- resmi[pos1]
  res[pos1na] <- resmi[pos1na] <- 0
  names(res)[pos1na] <- names(resmi)[pos1na] <- ed[pos1na]
  mi <- mi[pos2]
  mi[pos2na] <- 0
  names(mi)[pos2na] <- ed[pos2na]
  res <- res + sign(mi)
  resmi <- resmi + mi
}
resmi <- resmi / res

# 合并边数期望以 Poisson 上尾检验评估 bootstrap 稳健性并做多重校正
regs <- readLines("Data/regulators.txt")
genes <- sapply(strsplit(readLines("Data/expression_genes.txt")[-1], "\t"), function(x) x[1])
tedges <- length(genes) * length(which(genes %in% regs))
lambda <- sum(edges) / tedges
pvals <- p.adjust(ppois(res, lambda, lower.tail = FALSE), "fdr")
sig <- sort(pvals[pvals < 0.05])

net_tab <- cbind(
  t(sapply(strsplit(names(sig), collapse.key), function(x) x)),
  signif(resmi[match(names(sig), names(resmi))], 5),
  signif(sig, 5)
)
write.table(net_tab, "results/finalNetwork_4col.tsv", sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

# 网络转 Regulon 后以 VIPER 推断各样本调控子（蛋白质）活性
dset <- as.matrix(readRDS("Data/expression_matrix.rds"))
regulon <- aracne2regulon(afile = "results/finalNetwork_4col.tsv", eset = dset, format = "3col")
vp <- viper(dset, regulon, method = "none")

saveRDS(vp, "results/viper_protein_activity.rds")
