import urllib.request
import os

print("⏬ 1. 正在从 RCSB PDB 下载 EGFR C797S 蛋白结构 (6S9C)...")
url = "https://files.rcsb.org/download/6S9C.pdb"
urllib.request.urlretrieve(url, "6s9c_raw.pdb")

print("🧹 2. 正在清洗蛋白结构（提取 A 链蛋白，去除水分子及杂质）...")
with open("6s9c_raw.pdb", "r") as f_in, open("receptor_clean.pdb", "w") as f_out:
    for line in f_in:
        # 保留 A 链的标准蛋白残基原子 (ATOM)
        if line.startswith("ATOM") and line[21] == 'A':
            f_out.write(line)

print("✅ 受体蛋白清洗完成，生成 `receptor_clean.pdb`！")