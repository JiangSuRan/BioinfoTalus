## dbDiffusion VAE 潜空间 classifier-free 引导扩散采样、表达解码与真实扰动评估
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import pearsonr
import torch

sys.path.append("functions/diffusion/classifier_free_v3_no_saver")
from VAE.VAE_model import VAE
from guided_diffusion.script_util_classifier_free_v3 import (
    create_model_and_diffusion,
    model_and_diffusion_defaults,
)

device = torch.device("cuda")
torch.manual_seed(0)

adata = sc.read_h5ad("Data/perturb_lognorm.h5ad")
emb_dict = np.load("Data/emb_dict.npy", allow_pickle=True).item()
target_pert = "ZC3H13"
ctrl_label = "non-targeting"

vae = VAE(num_genes=adata.n_vars, device=str(device), hidden_dim=128)
vae.load_state_dict(torch.load("Data/model/vae.pt", map_location=device, weights_only=False))
vae.to(device).eval()

cfg = model_and_diffusion_defaults()
cfg["input_dim"] = 128
cfg["class_input_dim"] = len(emb_dict[target_pert])
cfg["eta"] = 1.5
cfg["all_mask"] = False

model, diffusion = create_model_and_diffusion(**cfg)
model.load_state_dict(torch.load("Data/model/trained_Unet_fullLoop800000_v3.pt", map_location=device, weights_only=False))
model.to(device).eval()

# classifier-free 采样内部复制条件拼接会限制单批容量，按 100 细胞分批累计采样潜向量后统一 VAE 解码
n_samples = 5000
safe_batch = 100
emb_vec = np.asarray(emb_dict[target_pert], dtype=np.float32)
latent_chunks = []
got = 0
while got < n_samples:
    b = min(safe_batch, n_samples - got)
    cond = torch.from_numpy(emb_vec).unsqueeze(0).repeat(b, 1).to(device)
    with torch.no_grad():
        sample, _ = diffusion.p_sample_loop(model, (b, 128), cond, clip_denoised=False, model_kwargs={})
    latent_chunks.append(sample.cpu().numpy())
    got += b

latent = np.concatenate(latent_chunks, axis=0)
decoded_chunks = []
for i in range(0, len(latent), 256):
    z = torch.tensor(latent[i:i + 256], dtype=torch.float32, device=device)
    with torch.no_grad():
        decoded_chunks.append(vae.forward(z, return_decoded=True).cpu().numpy())

gen_expr = np.concatenate(decoded_chunks, axis=0)
np.save("results/ZC3H13_generated.npy", gen_expr)

X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
labels = adata.obs["perturbation"].astype(str).to_numpy()
ctrl_mean = X[labels == ctrl_label].mean(axis=0)
real_mean = X[labels == target_pert].mean(axis=0)
gen_mean = gen_expr.mean(axis=0)

r_mean = pearsonr(real_mean, gen_mean)[0]
r_eff = pearsonr(real_mean - ctrl_mean, gen_mean - ctrl_mean)[0]
pd.DataFrame([{
    "target": target_pert,
    "n_generated": gen_expr.shape[0],
    "pearson_mean_expr": float(r_mean),
    "pearson_effect_size": float(r_eff),
}]).to_csv("results/dbdiffusion_eval_metrics.csv", index=False)
