
# 03_vector_indexing.py
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs

# 加载奥希替尼及四代 EGFR 候选分子 (如 BLU-945 骨架)
OSIMERTINIB_SMILES = "CN1CCN(CC1)c2cc(c(cc2Nc3nc(cn3)c4cn(c5ccccc45)C)NC(=O)C=C)OC"

def smiles_to_bitvector(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return fp

def calculate_tanimoto_similarity():
    ref_fp = smiles_to_bitvector(OSIMERTINIB_SMILES)
    df = pd.read_csv("egfr_filtered_druglike.csv").head(200) # 取前200个做高精度计算
    
    similarities = []
    for smiles in df['smiles']:
        fp = smiles_to_bitvector(smiles)
        if fp and ref_fp:
            score = DataStructs.TanimotoSimilarity(ref_fp, fp)
            similarities.append(score)
        else:
            similarities.append(0.0)
            
    df['tanimoto_similarity'] = similarities
    df = df.sort_values(by='tanimoto_similarity', ascending=False)
    df.to_csv("egfr_similarity_ranked.csv", index=False)
    print("✅ 基于 Tanimoto 相似度算法排序完成，顶尖候选分子已生成！")

if __name__ == "__main__":
    calculate_tanimoto_similarity()