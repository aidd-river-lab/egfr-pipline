import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def apply_lipinski_filter(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    mw = Descriptors.MolWt(mol)           # 分子量 ≤ 500
    logp = Descriptors.MolLogP(mol)        # 脂水分配系数 ≤ 5
    hbd = Descriptors.NumHDonors(mol)      # 氢键供体 ≤ 5
    hba = Descriptors.NumHAcceptors(mol)   # 氢键受体 ≤ 10
    
    is_druglike = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
    return is_druglike, mw, logp, hbd, hba

def process_dataset():
    print("🚀 [Step 2/3] 开始读取数据并执行 Lipinski Rule of 5 过滤...")
    df = pd.read_csv("egfr_chembl_compounds.csv")
    print(f"📥 待分析的数据总量: {len(df)} 条")
    
    results = []
    for idx, row in df.iterrows():
        filter_res = apply_lipinski_filter(row['canonical_smiles'])
        if filter_res and filter_res[0]:  # 仅保留符合类药性规范的分子
            results.append({
                'chembl_id': row['molecule_chembl_id'],
                'smiles': row['canonical_smiles'],
                'ic50_nm': row['standard_value'],
                'mw': round(filter_res[1], 2),
                'logp': round(filter_res[2], 2)
            })
            
    filtered_df = pd.DataFrame(results)
    output_csv = "egfr_filtered_druglike.csv"
    filtered_df.to_csv(output_csv, index=False)
    print(f"✅ 理化性质过滤完成！符合 Druglike 规范的候选分子: {len(filtered_df)} 个")
    print(f"💾 结果已保存至 `{output_csv}`\n")

if __name__ == "__main__":
    process_dataset()