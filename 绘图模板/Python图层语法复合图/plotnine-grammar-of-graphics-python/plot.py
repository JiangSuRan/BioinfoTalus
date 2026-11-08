## plotnine Python 图层语法分面统计图（ggplot 风格出版级输出）
import numpy as np
import pandas as pd
from plotnine import ggplot, aes, geom_bar, geom_hline, geom_text, theme_tufte, theme, coord_flip, scale_fill_manual, ggsize, ggsave, element_text, options

adata_obs = pd.read_csv("Data/peak_obs.csv")

# 比例条形图：按注释分组的填充占比 + 参考线 + 标签层
label_df = pd.DataFrame({
    "cluster": adata_obs["final annotation"].unique(),
    "lab": "n",
})

p = (
    ggplot(adata_obs, aes(x="final annotation", fill="type"))
    + geom_bar(position="fill")
    + theme_tufte()
    + theme(dpi=200, figure_size=(3, 8), axis_text_x=element_text(rotation=90))
    + scale_fill_manual(np.unique(adata_obs["type"].astype(str)).tolist())
    + coord_flip()
    + geom_hline(
        yintercept=np.sum(adata_obs["type"].astype(str).str.endswith("_cl13")) / adata_obs.shape[0],
        linetype="dotted",
        color="gray",
    )
    + geom_text(aes(label="lab", x="cluster"), y=1.02, inherit_aes=False, data=label_df, size=8)
    + ggsize(300, 800)
)

ggsave(p, "results/plotnine_barplot.pdf", verbose=False)
