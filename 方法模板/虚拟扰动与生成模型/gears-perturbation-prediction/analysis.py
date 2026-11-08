## GEARS 图神经网络单基因与组合基因遗传扰动响应预测
import pickle

from gears import PertData, GEARS

pert_data = PertData("./Data/gears_data")
pert_data.load(data_name="norman")
pert_data.prepare_split(split="simulation", seed=1)
pert_data.get_dataloader(batch_size=32, test_batch_size=128)

gears_model = GEARS(pert_data, device="cuda", weight_bias_track=False,
                    proj_name="pertnet", exp_name="pertnet")
gears_model.model_initialize(hidden_size=64)
gears_model.train(epochs=20, lr=1e-3)
gears_model.save_model("results/gears_model")

# 预测未见过的单基因与双基因组合扰动的转录响应（含差异基因与表型得分）
pred = gears_model.predict([["FEV"], ["FEV", "AHR"]])
with open("results/gears_predictions.pkl", "wb") as f:
    pickle.dump(pred, f)
