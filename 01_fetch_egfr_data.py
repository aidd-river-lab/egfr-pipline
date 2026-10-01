import pandas as pd
from chembl_webresource_client.new_client import new_client

def fetch_egfr_compounds():
    print("🚀 1. 正在连接 ChEMBL API 检索 Target: EGFR (CHEMBL203)...")
    activity = new_client.activity
    
    # 过滤单位为 nM 且测量指标为 IC50 的记录
    # 加 limit(1500) 避免网络卡死，快速获取前 1500 条高质量数据
    res = activity.filter(
        target_chembl_id="CHEMBL203",
        standard_type="IC50",
        standard_units="nM"
    ).only(['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'pchembl_value'])[:1500]
    
    print("⏬ 2. 正在快速下载 1500 条活性数据并解析...")
    df = pd.DataFrame(res)
    
    print("🧹 3. 正在清洗与过滤空值...")
    df = df.dropna(subset=['canonical_smiles', 'standard_value'])
    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    
    # 仅保留具有真实数值的小分子，按 IC50 活性（越小活性越高）排序
    df = df.dropna(subset=['standard_value'])
    df = df.sort_values(by='standard_value', ascending=True)
    
    output_csv = "egfr_chembl_compounds.csv"
    df.to_csv(output_csv, index=False)
    print(f"✅ 数据提取完成！共快速获取 {len(df)} 条有效数据，已保存至 `{output_csv}`。")

if __name__ == "__main__":
    fetch_egfr_compounds()