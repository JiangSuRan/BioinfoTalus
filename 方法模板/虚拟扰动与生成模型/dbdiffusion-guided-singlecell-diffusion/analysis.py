## dbDiffusion VAE 潜空间 classifier-free 引导扩散单细胞生成采样
import os

import numpy as np
import torch as th

import sys
sys.path.append("functions/diffusion/classifier_free_v3_no_saver")
from guided_diffusion import dist_util, logger
from guided_diffusion.script_util_classifier_free_v3_acc import (
    create_model_and_diffusion,
    model_and_diffusion_defaults,
    args_to_dict,
)

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
device = th.device("cuda" if th.cuda.is_available() else "cpu")

# 固定随机种子
random_seed = 1234
np.random.seed(random_seed)
th.manual_seed(random_seed)

# 类别（细胞态/亚群）嵌入：VAE 编码后的潜向量，复制为批次维度
class_emb_for_one_class = np.load("Data/class_emb_target.npy")
pert_name_for_one_class = "target_state"
folder_path = "results/"
neighbor = "ref"

args_defaults = model_and_diffusion_defaults()
args_defaults.update(
    {
        "num_samples": 30000,
        "batch_size": 10000,
        "eta": 1.5,
        "all_mask": False,
        "model_path": "Data/model/trained_Unet_fullLoop800000_v3.pt",
        "sample_dir": f"{folder_path}/{neighbor}_v3_eta1.5_30000gen.npz",
    }
)
args = args_defaults

class_emb = th.from_numpy(class_emb_for_one_class)
class_emb_to_use = class_emb.unsqueeze(0).repeat(args["batch_size"], 1).to(th.float32).to(dist_util.dev())
args["class_input_dim"] = class_emb_to_use.shape[1]

dist_util.setup_dist()
logger.configure(dir="results/sample_logs")

model, diffusion = create_model_and_diffusion(**args_to_dict(args, model_and_diffusion_defaults().keys()))
model.load_state_dict(dist_util.load_state_dict(args["model_path"], map_location="cpu"))
model.to(dist_util.dev())
model.eval()

# DDIM 采样（加速版）：以类别嵌入条件生成 VAE 潜空间样本
all_cells = []
while len(all_cells) * args["batch_size"] < args["num_samples"]:
    sample_fn = diffusion.p_sample_loop if not args["use_ddim"] else diffusion.ddim_sample_loop
    sample, traj = sample_fn(
        model,
        (args["batch_size"], args["input_dim"]),
        class_emb_to_use,
        clip_denoised=args["clip_denoised"],
        model_kwargs={},
        start_time=diffusion.betas.shape[0],
        accel_flag=True,
    )
    all_cells.append(sample.cpu().numpy())

arr = np.concatenate(all_cells, axis=0)
np.savez(args["sample_dir"], cell_gen=arr)
