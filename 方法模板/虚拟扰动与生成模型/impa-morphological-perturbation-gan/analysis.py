## IMPA 风格迁移生成模型预测药物扰动的细胞形态图像表型
import yaml
import torch
from omegaconf import OmegaConf

from IMPA.dataset.data_loader import CellDataLoader
from IMPA.solver import IMPAmodule

path_to_config = "Data/impa_config/bbbc021_six.yaml"

# 读取并修正相对路径后装配配置
with open(path_to_config, "r") as file:
    config = yaml.safe_load(file)

config["image_path"] = "Data/" + config["image_path"]
config["data_index_path"] = "Data/" + config["data_index_path"]
config["embedding_path"] = "Data/" + config["embedding_path"]
args = OmegaConf.create(config)

dataloader = CellDataLoader(args)

# 加载训练检查点（200 epoch）
checkpoint_dir = "Data/checkpoints/bbbc021_six/"
solver = IMPAmodule(args, checkpoint_dir, dataloader)
solver._load_checkpoint(200)

val_dataloader = dataloader.val_dataloader()
mol2id = dataloader.mol2id

# 对未见药物（嵌入向量）执行图像翻译并保存变换后的细胞形态表征
from tqdm import tqdm

transformed = {}
controls = []
for batch in tqdm(val_dataloader):
    _, _, _, y_hat = solver.forward(batch)
    for i, mol_id in enumerate(batch["mol_onehot"].argmax(dim=1).tolist()):
        mol_name = [k for k, v in mol2id.items() if v == mol_id][0]
        transformed.setdefault(mol_name, []).append(y_hat[i].detach().cpu())

torch.save(transformed, "results/impa_transformed_phenotypes.pt")
