## cNMF 基因表达程序（GEP）汇总与已知共表达模块 Jaccard 相似性比对
fns <- list.files(pattern = "*.csv", path = "Data/GEPs")
geplist <- list()
for (fn in fns) {
  # 文件名前缀标识细胞类；CSV 首列为基因名，其余每列为一个 GEP 的基因列表
  file.pref <- sub("_.*", "", fn)
  df <- read.delim(file.path("Data/GEPs", fn), sep = ",")
  df <- df[, -1]
  ll <- lapply(colnames(df), function(x) c(df[, x]))
  names(ll) <- colnames(df)
  geplist[[file.pref]] <- ll
}

saveRDS(geplist, "results/collated_GEPlists.rds")

coremap.list <- readRDS("Data/modules_resource_coremap.rds")

# pref 为 GEP 文件前缀（细胞类），ccpref 为参考模块映射的键；ccpref = "All" 时汇总全部细胞类的参考模块
pref <- "Ast"
ccpref <- "Ast"
ll <- geplist[[pref]]
usemap.list <- coremap.list

# 参考映射为 gene -> module 编号的命名向量，先翻转成 module -> genes 列表
if (ccpref != "All") {
  usemap <- usemap.list[[ccpref]]
  modnums <- sort(unique(usemap))
  cc <- lapply(modnums, function(x) names(usemap)[usemap == x])
  names(cc) <- paste0("M", modnums)
} else {
  cc <- list()
  for (name in names(usemap.list)) {
    usemap <- usemap.list[[name]]
    modnums <- sort(unique(usemap))
    subcc <- lapply(modnums, function(x) names(usemap)[usemap == x])
    names(subcc) <- paste0(name, "-M", modnums)
    cc <- c(cc, subcc)
  }
}

NL <- length(ll)
NC <- length(cc)
imat <- matrix(0, nrow = NL, ncol = NC)
for (i in 1:NL) {
  gep <- ll[[i]]
  for (j in 1:NC) {
    imat[i, j] <- length(intersect(gep, cc[[j]]))
  }
}

# Jaccard = 交集 / 并集，并集大小由两集合大小之和减交集得到
lct <- sapply(ll, length)
cct <- sapply(cc, length)
umat <- (matrix(rep(cct, NL), byrow = TRUE, nrow = NL, ncol = NC) +
  matrix(rep(lct, NC), byrow = FALSE, nrow = NL, ncol = NC)) - imat
jmat <- imat / umat
colnames(jmat) <- names(cc)
rownames(jmat) <- names(ll)
colnames(imat) <- names(cc)
rownames(imat) <- names(ll)

saveRDS(list(intersect = imat, jaccard = jmat), "results/gep_module_similarity.rds")
