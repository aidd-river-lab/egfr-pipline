import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs

# 1. 参照分子：奥希替尼 (Osimertinib) 标准 SMILES
OSIMERTINIB_SMILES = "CN1CCN(CC1)c2cc(c(cc2Nc3nc(cn3)c4cn(c5ccccc45)C)NC(=O)C=C)OC"

def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    # 提取 2048 位 Morgan 分子指纹 (Radius = 2)
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

def run_similarity_search():
    print("🚀 [Step 3/4] 正在计算与奥希替尼的 Tanimoto 结构相似度...")
    ref_fp = smiles_to_fp(OSIMERTINIB_SMILES)
    
    df = pd.read_csv("egfr_knowledge_filtered.csv")
    
    similarities = []
    for smiles in df['smiles']:
        fp = smiles_to_fp(smiles)
        if fp and ref_fp:
            score = DataStructs.TanimotoSimilarity(ref_fp, fp)
            similarities.append(round(score, 4))
        else:
            similarities.append(0.0)
            
    df['tanimoto_similarity'] = similarities
    # 按相似度由高到低排序
    df = df.sort_values(by='tanimoto_similarity', ascending=False)
    
    # 选取最相似的前 5 个候选分子，保存供下一步对接
    top_candidates = df.head(5)
    output_csv = "egfr_top_candidates.csv"
    top_candidates.to_csv(output_csv, index=False)
    
    print(f"✅ 筛选完成！已选出 Top 5 最相似的候选分子，结果已保存至 `{output_csv}`：\n")
    print(top_candidates[['chembl_id', 'ic50_nm', 'mw', 'tanimoto_similarity']])

if __name__ == "__main__":
    run_similarity_search()