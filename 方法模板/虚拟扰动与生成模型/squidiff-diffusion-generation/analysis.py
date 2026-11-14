## Squidiff 扩散模型单细胞转录状态条件生成与扰动采样
import numpy as np
import torch

from Squidiff import dist_util
from Squidiff.script_util import (
    create_model_and_diffusion,
    model_and_diffusion_defaults,
)

# 按训练配置重建模型与扩散过程并载入权重
default_args = model_and_diffusion_defaults()
default_args.update(
    {
        "schedule_sampler": "uniform",
        "batch_size": 64,
        "use_fp16": True,
        "gene_size": 5000,
        "output_dim": 512,
        "num_layers": 3,
        "class_cond": False,
        "use_encoder": True,
        "use_ddim": True,
        "diffusion_steps": 1000,
        "model_path": "Data/squidiff_model.pt",
        "use_drug_structure": True,
        "comb_num": 1,
        "drug_dimension": 1024,
    }
)
model, diffusion = create_model_and_diffusion(**default_args)
model.load_state_dict(dist_util.load_state_dict(default_args["model_path"]))
model.to(dist_util.dev())
model.eval()

# 条件向量 z_sem（如药物结构编码均值）经 DDIM 反向采样生成细胞转录状态
z_sem = torch.tensor(np.load("Data/condition_latent.npy"), dtype=torch.float32, device=dist_util.dev())
sample_fn = diffusion.ddim_sample_loop if default_args["use_ddim"] else diffusion.p_sample_loop
samples = sample_fn(
    model,
    shape=(z_sem.shape[0], default_args["gene_size"]),
    model_kwargs={"z_mod": z_sem},
    noise=None,
)

np.save("results/squidiff_samples.npy", samples.detach().cpu().numpy())
