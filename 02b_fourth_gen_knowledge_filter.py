# 02b_fourth_gen_knowledge_filter.py
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog

# ===========================================================================
# 1. 踩坑黑名单字典 (C797S 特定毒性/代谢警示结构)
# ===========================================================================
PITFALL_SMARTS = {
    # 坑 1: BLU-945 式未掩蔽苯胺 (极易代谢生成醌亚胺，引发肝毒性)
    "blu945_unmasked_aniline": "c1ccccc1[NH2,NH1CH3]",
    
    # 坑 2: 裸露呋喃/噻吩 (极易氧化开环生成毒性二醛/S-氧化物)
    "reactive_furan_thiophene": "o1cccc1,s1cccc1",
    
    # 坑 3: BBT-176 式 hERG 心脏毒性核心药效团 (强碱性长链中心)
    "herg_basic_amine_chain": "CCN(CC)CCN",
    
    # 坑 4: 易发生脱甲基/快速脱毒代谢的软位点
    "labile_methoxy": "c1ccccc1OC",
}

# 加载 RDKit 内置的通用 PAINS 假阳性/毒性警示库
pains_params = FilterCatalog.FilterCatalogParams()
pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
pains_catalog = FilterCatalog.FilterCatalogCatalog(pains_params)

# ===========================================================================
# 2. 优势白名单字典 (四代药抓牢 Ser797 丝氨酸的关键氢键锚点)
# ===========================================================================
SER797_ANCHOR_SMARTS = [
    "S(=O)(=O)",        # 磺酰基 (强氢键受体)
    "c1cnccc1",         # 吡啶环 (氮原子受体)
    "n1c(N)nc(N)cc1"    # 类似四代分子的氢键网骨架
]

def check_toxic_and_pitfalls(mol):
    """黑名单检测：包含结构坑、PAINS 和理化性质坑"""
    hits = []
    
    # A. 检查自定义四代踩坑结构
    for name, smarts_str in PITFALL_SMARTS.items():
        for smarts in smarts_str.split(","):
            patt = Chem.MolFromSmarts(smarts)
            if patt and mol.HasSubstructMatch(patt):
                hits.append(f"命中结构坑: {name}")
                break
                
    # B. 检查 PAINS 官方库
    if pains_catalog.HasMatch(mol):
        hits.append("命中 PAINS 假阳性库")
        
    # C. 理化性质踩坑检测 (结合 BLU / BBT 临床教训)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if logp > 3.8 and mw > 500:
        hits.append("理化坑: 高脂溶性+大分子量 (肝脏代谢极快/易积聚肝毒性)")
    if logp > 3.5 and tpsa < 70 and rotatable_bonds > 6:
        hits.append("心脏毒性坑: 高 hERG 阻断风险")
        
    return hits

def calculate_ser797_anchor_score(mol):
    """白名单检测：计算 Ser797 锚定得分"""
    score = 0
    for smarts in SER797_ANCHOR_SMARTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            score += 1
    return score

def run_knowledge_pipeline():
    print("🚀 启动 [四代 EGFR-C797S 知识库融合过滤引擎]...")
    
    input_csv = "egfr_filtered_druglike.csv"
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"❌ 未找到 `{input_csv}`，请先运行 02_cheminformatics_filter.py！")
        return

    results = []
    for idx, row in df.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
            
        pitfalls = check_toxic_and_pitfalls(mol)
        
        # 硬性拦截：踩中结构坑或重度理化坑的直接淘汰
        if len(pitfalls) > 0:
            continue
            
        anchor_score = calculate_ser797_anchor_score(mol)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        results.append({
            'chembl_id': row['chembl_id'],
            'smiles': smiles,
            'ic50_nm': row['ic50_nm'],
            'mw': round(mw, 2),
            'logp': round(logp, 2),
            'ser797_anchor_score': anchor_score
        })
        
    res_df = pd.DataFrame(results)
    
    # 优先保留具备更多 Ser797 锚定特征的优质分子
    res_df = res_df.sort_values(by=['ser797_anchor_score', 'ic50_nm'], ascending=[False, True])
    
    output_csv = "egfr_knowledge_filtered.csv"
    res_df.to_csv(output_csv, index=False)
    
    print(f"✅ 知识库过滤完成！成功规避历史毒性坑，保留了 {len(res_df)} 个高质量候选分子。")
    print(f"💾 结果已保存至 `{output_csv}`\n")

if __name__ == "__main__":
    run_knowledge_pipeline()