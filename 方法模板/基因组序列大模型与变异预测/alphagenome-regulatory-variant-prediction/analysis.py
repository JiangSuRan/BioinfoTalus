## AlphaGenome 基因组序列调控效应预测（RNA-seq/ATAC/CHIP/接触图与变异打分）
import os

import numpy as np
import pandas as pd

os.environ["XLA_FLAGS"] = " ".join([
    "--xla_gpu_deterministic_ops",
    "--xla_gpu_enable_scatter_determinism_expander=True",
    "--xla_gpu_enable_triton_gemm=False",
])
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.9"

from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome_research.model import dna_model

# 加载多折模型（Kaggle 下载）与人类 track 元数据
model = dna_model.create_from_kaggle("all_folds")
metadata = model.output_metadata(dna_model.Organism.HOMO_SAPIENS)

# GTF 转 MANE Select 转录本提取器（可视化注释）
gtf = pd.read_feather(
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)
gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
gtf_transcripts = gene_annotation.filter_to_mane_select_transcript(gtf_transcripts)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)

# 区间预测：1Mb 窗口的多组学 track（RNA-seq / ATAC / CHIP / 接触图）
interval = genome.Interval("chr19", 1_058_7331, 1_163_5907, strand="+")
predictions = model.predict_interval(
    interval,
    requested_outputs={
        dna_model.OutputType.RNA_SEQ,
        dna_model.OutputType.DNASE,
        dna_model.OutputType.ATAC,
        dna_model.OutputType.CHIP_TF,
        dna_model.OutputType.CHIP_HISTONE,
        dna_model.OutputType.CONTACT_MAPS,
    },
    ontology_terms=["EFO:0001187"],
)

# 变异效应打分：参考/替代序列 1Mb 窗口预测差
variant = genome.Variant.from_str("chr21:46126238:G>C")
variant_interval = variant.reference_interval.resize(2**20)
variant_predictions = model.predict_variant(
    interval=variant_interval,
    variant=variant,
    requested_outputs={dna_model.OutputType.RNA_SEQ},
    ontology_terms=["EFO:0001187"],
)

# 汇总 RNA-seq 变异效应到 csv
ref_alt = variant_predictions.rna_seq.values
variant_effect = pd.DataFrame(
    {"track": np.arange(ref_alt.shape[-1]), "delta": ref_alt.reshape(ref_alt.shape[0], -1).ravel()}
)
variant_effect.to_csv("results/alphagenome_variant_effect.csv", index=False)
