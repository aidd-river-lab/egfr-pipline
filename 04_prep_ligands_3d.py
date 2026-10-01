import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_3d_structures():
    print("🚀 [Step 4/4] 开始将候选分子的 2D SMILES 转化为 3D 构象文件 (.pdb)...")
    
    # 创建存放 3D 结构的文件夹
    os.makedirs("ligands_3d", exist_ok=True)
    
    df = pd.read_csv("egfr_top_candidates.csv")
    
    for idx, row in df.iterrows():
        chembl_id = row['chembl_id']
        smiles = row['smiles']
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
            
        # 1. 加氢（显式还原氢原子位置）
        mol = Chem.AddHs(mol)
        
        # 2. 使用 ETKDGv3 算法生成 3D 空间坐标
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        embed_status = AllChem.EmbedMolecule(mol, params)
        
        if embed_status == 0: # 生成成功
            # 3. MMFF94 力场能量最小化（消除原子碰撞，优化构象）
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except:
                pass
                
            # 4. 导出为 PDB 文件
            pdb_path = f"ligands_3d/{chembl_id}.pdb"
            Chem.MolToPDBFile(mol, pdb_path)
            print(f"  └─ ✅ 成功生成 {chembl_id} 的 3D 构象 ➔ `{pdb_path}`")
        else:
            print(f"  └─ ❌ {chembl_id} 3D 构象生成失败")

    print("\n🎉 所有候选分子的 3D 结构已就绪，位于 `ligands_3d/` 目录下！")

if __name__ == "__main__":
    generate_3d_structures()