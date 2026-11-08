## coolbox 多组学基因组浏览器轨迹图（BigWig + 基因结构 + 高亮区域）
import numpy as np
import pandas as pd
from coolbox.api import *

# 目标 peak 区域解析为 chr:start-end
peak_info = pd.read_csv("Data/peaks.more_info.bed", sep="\t", header=None)
new_regions = {"T cell activation arch": "peak.44945"}
def region_str(df):
    return "%s:%s-%s" % (df[0].values[0], df[1].values[0], df[2].values[0])
loc = region_str(peak_info[peak_info[3] == new_regions["T cell activation arch"]])

# BigWig 轨迹（各细胞亚群假批次 track）+ 基因结构 + peak atlas
bws = [
    ("Data/bigwig/Cd4_Tfh.bw", "#bd6960"),
    ("Data/bigwig/Cd8_dys.bw", "#0039d6"),
    ("Data/bigwig/Cd8_prog.bw", "#abbadf"),
]

frame = None
for bw_path, color in bws:
    track = BigWig(bw_path) + Title(bw_path.split("/")[-1].split(".")[0]) + Color(color)
    frame = track if frame is None else frame + track

frame = (
    frame
    + XAxis()
    + BED("Data/annotation/sorted.gencode.vM25.bed", gene_style="flybase")
    + BED("Data/annotation/peak_atlas.bed")
)

# 渲染并保存矢量图（Browser.save 支持指定范围）
bsr = Browser(frame, init_range=loc)
bsr.save("results/coolbox_multitrack.pdf", dpi=300)
