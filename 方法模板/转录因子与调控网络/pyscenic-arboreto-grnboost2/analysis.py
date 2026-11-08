## pySCENIC 三步管线（grnboost2 共表达 + cisTarget 修剪 + AUCell 活性）
import glob
import os
import pickle

import scanpy as sc
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.aucell import aucell
from pyscenic.prune import df2regulons, prune2df
from pyscenic.utils import modules_from_adjacencies

DATA_FOLDER = "results/scenic"
RESOURCES_FOLDER = "Data/scenic"
DATABASES_GLOB = os.path.join(RESOURCES_FOLDER, "*feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
HG_TFS_FNAME = os.path.join(RESOURCES_FOLDER, "allTFs_hg38.txt")
MIN_GENES = 20

os.makedirs(DATA_FOLDER, exist_ok=True)
tf_names = load_tf_names(HG_TFS_FNAME)

db_fnames = glob.glob(DATABASES_GLOB)
dbs = [RankingDatabase(fname=fname, name=os.path.splitext(os.path.basename(fname))[0]) for fname in db_fnames]

adata = sc.read_h5ad("Data/myeloid_counts.h5ad")
ex_matrix = adata.to_df()

# Phase 1：grnboost2 共表达邻接 + 共表达模块
adjacencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=False)
adjacencies.to_csv(f"{DATA_FOLDER}/adjacencies.tsv", index=False, sep="\t")
modules = list(modules_from_adjacencies(adjacencies, ex_matrix, min_genes=MIN_GENES))
with open(f"{DATA_FOLDER}/modules.p", "wb") as f:
    pickle.dump(modules, f)

# Phase 2：cisTarget motif 修剪为 regulons
df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, num_workers=6)
df.to_csv(f"{DATA_FOLDER}/motifs.csv")
regulons = df2regulons(df)
with open(f"{DATA_FOLDER}/regulons.p", "wb") as f:
    pickle.dump(regulons, f)

# Phase 3：AUCell regulon 活性矩阵
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
auc_mtx.to_csv(f"{DATA_FOLDER}/aucell.csv")

sc.settings.verbosity = 0
adata.obs = adata.obs.join(auc_mtx, how="left")
adata.write_h5ad("results/pyscenic_annotated.h5ad")
