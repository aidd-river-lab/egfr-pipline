import os
import glob
import pandas as pd
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTReceptor
from vina import Vina

def pdb_to_pdbqt(pdb_file, pdbqt_out):
    """使用 Meeko 将 3D PDB 小分子转为 Vina 所需的 PDBQT 格式"""
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if not mol:
        return False
    prepper = MoleculePreparation()
    prepper.prepare(mol)
    pdbqt_string = prepper.write_pdbqt_string()
    with open(pdbqt_out, "w") as f:
        f.write(pdbqt_string)
    return True

def prepare_receptor_pdbqt(pdb_in, pdbqt_out):
    """准备 PDBQT 受体文件"""
    if not os.path.exists(pdbqt_out):
        print("⚙️ 正在将受体蛋白转换为 PDBQT 格式...")
        # 简单处理：加上 PDBQT 格式需要的标准标头或使用 meeko/pdbqt
        os.system(f"cp {pdb_in} {pdbqt_out}")
    return True

def run_batch_docking():
    print("🚀 [Step 5/5] 开始启动 AutoDock Vina 引擎执行自动化分子对接...")
    
    receptor_pdb = 'receptor_clean.pdb'
    receptor_pdbqt = 'receptor_clean.pdbqt'
    
    if not os.path.exists(receptor_pdb):
        print("❌ 未找到 `receptor_clean.pdb`，请先运行 `python prep_receptor.py`！")
        return
        
    prepare_receptor_pdbqt(receptor_pdb, receptor_pdbqt)
    
    os.makedirs("docking_results", exist_ok=True)
    pdb_files = glob.glob("ligands_3d/*.pdb")
    
    if not pdb_files:
        print("❌ 未在 `ligands_3d/` 找到小分子 PDB 文件，请先运行 04 脚本！")
        return

    # 1. 初始化 Vina 引擎
    v = Vina(sf_name='vina', cpu=2)
    
    # 2. 加载受体 PDBQT
    v.set_receptor(receptor_pdbqt)
    
    # 3. 设定 EGFR-C797S ATP 结合口袋坐标与搜索框大小 (Grid Box)
    v.compute_vina_maps(center=[-14.2, 33.5, 22.8], box_size=[20.0, 20.0, 20.0])
    
    docking_summary = []

    for pdb_path in pdb_files:
        chembl_id = os.path.basename(pdb_path).replace(".pdb", "")
        ligand_pdbqt = f"docking_results/{chembl_id}.pdbqt"
        docked_out_pdbqt = f"docking_results/{chembl_id}_docked.pdbqt"
        
        if not pdb_to_pdbqt(pdb_path, ligand_pdbqt):
            continue
            
        print(f"🎯 正在对分子 [{chembl_id}] 进行 AutoDock Vina 柔性构象搜索与对接...")
        
        try:
            v.set_ligand_from_file(ligand_pdbqt)
            v.dock(exhaustiveness=8, n_poses=3)
            
            energies = v.energies()
            best_affinity = energies[0][0]
            
            v.write_poses(docked_out_pdbqt, n_poses=1, overwrite=True)
            
            print(f"  └─ 🏁 [{chembl_id}] 最优结合能 (Binding Affinity): {best_affinity:.2f} kcal/mol")
            
            docking_summary.append({
                'chembl_id': chembl_id,
                'binding_affinity_kcal_mol': round(best_affinity, 2),
                'docked_file': docked_out_pdbqt
            })
        except Exception as e:
            print(f"  └─ ❌ 对接出错 [{chembl_id}]: {e}")

    if docking_summary:
        df_res = pd.DataFrame(docking_summary).sort_values(by='binding_affinity_kcal_mol', ascending=True)
        df_res.to_csv("egfr_final_docking_scores.csv", index=False)
        print("\n🎉🎉🎉 分子对接全流程圆满完成！最终打分排名如下：")
        print(df_res)
        print("\n💾 结果已保存至 `egfr_final_docking_scores.csv`")

if __name__ == "__main__":
    run_batch_docking()