## PULSAR 供体级多尺度表征提取与疾病 logistic 分类
import numpy as np
import scanpy as sc
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import train_test_split

from pulsar.model import PULSAR
from pulsar.utils import extract_donor_embeddings_from_h5ad

np.random.seed(42)
torch.manual_seed(42)

adata = sc.read_h5ad("Data/lupus_uce_adata.h5ad")

# 预训练 PULSAR（HuggingFace）bfloat16 推理
pulsar_model = PULSAR.from_pretrained("KuanP/pulsar-pbmc")
pulsar_model.eval()
pulsar_model.cuda()
pulsar_model = pulsar_model.to(torch.bfloat16)

# 每供体采样 1024 细胞编码为供体级嵌入
donor_embedding_collection = extract_donor_embeddings_from_h5ad(
    adata,
    donor_id_key="donor_id",
    model=pulsar_model,
    sample_cell_num=1024,
)

# 供体表型对齐与特征矩阵组装
donor_id_to_disease = adata.obs.groupby("donor_id", observed=True)["disease"].first().to_dict()

X = []
y = []
donor_ids = []
for donor_id, donor_data in donor_embedding_collection.items():
    X.append(donor_data["embedding"][0])
    y.append(donor_id_to_disease[donor_id])
    donor_ids.append(donor_id)
X = np.vstack(X)
y = np.array(y)
donor_ids = np.array(donor_ids)

# 分层拆分 + 平衡权重 logistic 分类
X_train, X_test, y_train, y_test, donors_train, donors_test = train_test_split(
    X, y, donor_ids, test_size=0.2, random_state=0, stratify=y
)
clf = LogisticRegression(max_iter=1000, random_state=42, class_weight="balanced")
clf.fit(X_train, y_train)

y_test_pred = clf.predict(X_test)
test_accuracy = accuracy_score(y_test, y_test_pred)

report = classification_report(y_test, y_test_pred, zero_division=0)
with open("results/pulsar_classification_report.txt", "w") as f:
    f.write(f"Test accuracy: {test_accuracy:.4f}\n\n{report}")
np.savez("results/pulsar_donor_embeddings.npz", X=X, y=y, donor_ids=donor_ids)
