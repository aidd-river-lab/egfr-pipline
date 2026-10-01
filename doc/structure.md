# EGFR 药物虚拟筛选流程（egfr-pipline）项目说明

## 1. 项目是做什么的

这是一个**计算机辅助药物发现（CADD）的轻量级虚拟筛选流程**，目标靶点是 **EGFR（表皮生长因子受体）的 C797S 耐药突变体**。

背景：EGFR 是非小细胞肺癌（NSCLC）里最常见的驱动基因之一。第三代 EGFR 抑制剂奥希替尼（Osimertinib）上市后，临床上逐渐出现 **C797S 突变**导致的耐药——这个突变恰好发生在奥希替尼共价结合的半胱氨酸位点上，使药物失效。本项目就是围绕"能不能从已知对 EGFR 有活性的化合物里，筛出一批对 C797S 突变体仍然有潜力的候选分子"这个问题搭建的一条自动化流程。

整体思路是经典的虚拟筛选四步法：

```
已知活性数据(ChEMBL) → 类药性过滤 → 结构相似度排序 → 3D构象生成 → 分子对接打分
```

## 2. 整体流程图

```
prep_receptor.py                          01_fetch_egfr_data.py
(准备受体蛋白 6S9C)                         (从 ChEMBL 拉取 EGFR 活性数据)
        │                                           │
        │                                           ▼
        │                          02_cheminformatics_filter.py
        │                          (Lipinski 类药性五规则过滤)
        │                                           │
        │                                           ▼
        │                           03_similarity_search.py
        │                      (与奥希替尼做 Tanimoto 结构相似度排序，取 Top5)
        │                                           │
        │                                           ▼
        │                           04_prep_ligands_3d.py
        │                        (2D SMILES → 3D 构象，生成 .pdb)
        │                                           │
        └───────────────────┬───────────────────────┘
                             ▼
                   05_auto_docking.py
          (AutoDock Vina 批量对接，输出结合能排名)
```

`03_vector_indexing.py` **不在这条主流程里**，是一份内容和 `03_similarity_search.py` 高度相似的早期/实验版本（见下文第 5 节说明），实际产出链路用的是 `03_similarity_search.py`。

## 3. 各脚本详解

### 3.0 `prep_receptor.py` —— 准备受体蛋白

- **做什么**：从 RCSB PDB 数据库下载 EGFR C797S 突变体的晶体结构 **6S9C**，然后只保留 A 链的标准蛋白原子（`ATOM` 记录），去掉水分子、配体等杂质。
- **输入**：无（直接联网下载）。
- **输出**：
  - `6s9c_raw.pdb`：原始下载的完整 PDB 文件。
  - `receptor_clean.pdb`：清洗后只含 A 链蛋白原子的结构，供后续对接使用。
- **关键逻辑**：逐行过滤 `line.startswith("ATOM") and line[21] == 'A'`（PDB 格式第 22 列固定是链标识符）。

### 3.1 `01_fetch_egfr_data.py` —— 拉取 EGFR 生物活性数据

- **做什么**：通过 `chembl_webresource_client` 调用 ChEMBL 官方 API，检索靶点 **EGFR（CHEMBL203）** 下所有 **IC50（单位 nM）** 活性记录，取前 1500 条，清洗空值并按活性（IC50 越小越强）升序排序。
- **输出**：`egfr_chembl_compounds.csv`（字段：`canonical_smiles, molecule_chembl_id, pchembl_value, standard_value, value`）。
- **注意**：`limit(1500)` 是为了防止一次性拉取全部数据导致请求卡死，属于性能妥协，不是全量数据。

### 3.2 `02_cheminformatics_filter.py` —— 类药性过滤（Lipinski Rule of Five）

- **做什么**：用 RDKit 计算每个分子的四个理化指标，只保留同时满足以下条件的"类药"分子：
  - 分子量（MW）≤ 500
  - 脂水分配系数（LogP）≤ 5
  - 氢键供体数（HBD）≤ 5
  - 氢键受体数（HBA）≤ 10
- **输入**：`egfr_chembl_compounds.csv`
- **输出**：`egfr_filtered_druglike.csv`（字段：`chembl_id, smiles, ic50_nm, mw, logp`）
- 本次运行结果：1429 条输入 → 1025 条通过过滤。

### 3.3 `03_similarity_search.py` —— 结构相似度排序（主流程采用）

- **做什么**：以**奥希替尼（Osimertinib）**的 SMILES 为参照分子，计算每个候选分子的 **Morgan 指纹（半径 2，2048 位）**，用 **Tanimoto 系数**衡量结构相似度，按相似度从高到低排序，取 **Top 5** 作为下一步对接的候选分子。
- **输入**：`egfr_filtered_druglike.csv`
- **输出**：`egfr_top_candidates.csv`（在原字段基础上新增 `tanimoto_similarity`）
- **算法要点**：Morgan 指纹是一种基于原子邻域的圆形指纹（等价于 ECFP），Tanimoto 系数 = 两个指纹交集位数 / 并集位数，取值 0~1，越接近 1 说明结构越相似。

### 3.4 `03_vector_indexing.py` —— 相似度排序的早期/实验版本（非主流程）

- 和 3.3 几乎是同一套逻辑（同样的参照分子、同样的 Morgan 指纹 + Tanimoto），区别是：
  - 只取过滤后数据的**前 200 条**做计算（而不是全量）；
  - 不截取 Top 5，而是把全部排序结果写入 `egfr_similarity_ranked.csv`；
  - 没有被后续任何脚本读取（`04_prep_ligands_3d.py` 读的是 `egfr_top_candidates.csv`，不是 `egfr_similarity_ranked.csv`）。
- 建议：如果确认这是被主流程淘汰的草稿，可以删除或者移进一个 `archive/` 子目录，避免和 `03_similarity_search.py` 的编号冲突造成后来者误解"流程里有两个并行的第3步"。

### 3.5 `04_prep_ligands_3d.py` —— 2D 转 3D 构象

- **做什么**：把 Top 5 候选分子的 2D SMILES 转换成可用于分子对接的 3D 构象：
  1. `Chem.AddHs`：显式加氢；
  2. `AllChem.ETKDGv3`（`randomSeed=42`，结果可复现）：生成初始 3D 坐标；
  3. `AllChem.MMFFOptimizeMolecule`（MMFF94 力场，最多 500 次迭代）：能量最小化，消除原子碰撞；
  4. 导出为 `.pdb` 文件。
- **输入**：`egfr_top_candidates.csv`
- **输出**：`ligands_3d/<CHEMBL_ID>.pdb`（5 个文件，对应 Top 5 候选分子）

### 3.6 `05_auto_docking.py` —— AutoDock Vina 批量对接

- **做什么**：
  1. 把每个配体的 `.pdb` 用 **Meeko**（`MoleculePreparation`）转换成 Vina 需要的 `.pdbqt` 格式；
  2. 准备受体 `.pdbqt`（见下方"已知问题"，这一步目前是简化实现）；
  3. 用 **AutoDock Vina**（`vina` Python 绑定）设定对接口袋（`center=[-14.2, 33.5, 22.8]`，`box_size=[20, 20, 20]`，这组坐标是针对 6S9C 结构的 ATP 结合口袋硬编码的），对每个配体做柔性对接（`exhaustiveness=8`，输出 3 个构象取最优）；
  4. 汇总每个分子的最优结合能（kcal/mol，越负代表结合越强），按结合能升序排序。
- **输入**：`ligands_3d/*.pdb` + `receptor_clean.pdb`
- **输出**：
  - `docking_results/<CHEMBL_ID>.pdbqt`：配体 PDBQT
  - `docking_results/<CHEMBL_ID>_docked.pdbqt`：对接后构象
  - `egfr_final_docking_scores.csv`：最终排名结果（字段：`chembl_id, binding_affinity_kcal_mol, docked_file`）

## 4. 完整运行方式

按顺序执行（`prep_receptor.py` 和 `01_fetch_egfr_data.py` 互相独立，谁先跑都可以，但两者都要在 `05_auto_docking.py` 之前完成）：

```bash
python prep_receptor.py              # 准备受体蛋白
python 01_fetch_egfr_data.py         # 拉取 ChEMBL 活性数据
python 02_cheminformatics_filter.py  # Lipinski 过滤
python 03_similarity_search.py       # Tanimoto 相似度排序，取 Top5
python 04_prep_ligands_3d.py         # 生成 3D 构象
python 05_auto_docking.py            # 批量分子对接
```

## 5. 依赖环境

项目目前**没有 `requirements.txt`**，依赖需要根据脚本里的 `import` 手动安装。涉及到的第三方库：

| 库 | 用途 |
|---|---|
| `pandas` | 表格数据处理（CSV 读写） |
| `rdkit` | 分子解析、理化性质计算、指纹计算、3D 构象生成 |
| `chembl_webresource_client` | 访问 ChEMBL 官方 REST API |
| `numpy` | `03_vector_indexing.py` 用到 |
| `meeko` | 小分子 PDB → PDBQT 转换（为 Vina 对接做准备） |
| `vina` | AutoDock Vina 的 Python 绑定，执行分子对接 |

`urllib`、`os`、`glob` 为 Python 标准库，无需单独安装。

建议补一份 `requirements.txt`（参考同目录下 `cancer-egfr-landscape` 项目已有的写法）。

## 6. 产出文件一览

| 文件/目录 | 产出自 | 说明 |
|---|---|---|
| `6s9c_raw.pdb` | prep_receptor.py | 原始下载的 6S9C 晶体结构 |
| `receptor_clean.pdb` | prep_receptor.py | 清洗后只含 A 链蛋白原子 |
| `receptor_clean.pdbqt` | 05_auto_docking.py（运行时生成） | 受体 PDBQT（见已知问题） |
| `egfr_chembl_compounds.csv` | 01_fetch_egfr_data.py | ChEMBL 原始活性数据，1429 条 |
| `egfr_filtered_druglike.csv` | 02_cheminformatics_filter.py | 通过 Lipinski 过滤的分子，1025 条 |
| `egfr_top_candidates.csv` | 03_similarity_search.py | Top5 候选分子 + 相似度 |
| `egfr_similarity_ranked.csv` | 03_vector_indexing.py（非主流程） | 前 200 条的全量相似度排序 |
| `ligands_3d/*.pdb` | 04_prep_ligands_3d.py | Top5 候选分子的 3D 构象 |
| `docking_results/*.pdbqt` | 05_auto_docking.py | 配体 PDBQT + 对接后构象 |
| `egfr_final_docking_scores.csv` | 05_auto_docking.py | 最终对接打分排名 |

## 7. 已知问题 / 局限性

记录几个从当前代码和已产出的结果数据里能直接确认的问题，供后续修复参考：

1. **受体 PDBQT 转换目前是"假转换"，对接结果不可信。**
   `05_auto_docking.py` 里的 `prepare_receptor_pdbqt()` 实际只是 `cp receptor_clean.pdb receptor_clean.pdbqt`（已验证两个文件逐字节相同），并没有做真正的 PDBQT 转换（缺少 AutoDock 需要的原子类型、Gasteiger 部分电荷、可旋转键等信息，通常需要用 `prepare_receptor4.py`/`AutoDockTools` 或 Meeko 的受体准备工具完成）。
   对应地，仓库里现成的 `egfr_final_docking_scores.csv` 里 **5 个分子的结合能全部是 `0.0 kcal/mol`**——这不是一个真实的 Vina 打分结果（真实结合能通常是 -5 ~ -12 kcal/mol 量级的负数），说明这一步对接实际上没有产出有科学意义的结果，需要先修好受体准备这一步，再重新跑 05 脚本。

2. **当前 Top5 候选分子与奥希替尼的 Tanimoto 相似度全部是 0.0。**
   `egfr_top_candidates.csv` 里 5 行的 `tanimoto_similarity` 全是 `0.0`，意味着这批数据里没有任何分子在 Morgan 指纹空间下与奥希替尼有可检测的结构相似性。这有可能是真实现象（这批 ChEMBL 检索到的高活性分子骨架本来就和奥希替尼差异很大），但考虑到是"全部恰好为 0"而不是"普遍偏低"，建议之后跑一次时额外抽查几个中间结果，确认指纹计算/比较逻辑没有问题。

3. **`03_vector_indexing.py` 和 `03_similarity_search.py` 编号重复、职责重叠**，容易让人误以为是并行的两步。建议确认是否为迭代过程中留下的草稿，是的话归档或删除。

4. **没有 `requirements.txt`**，环境复现依赖阅读代码里的 `import` 手动安装，门槛较高（尤其是 `vina`/`meeko` 这类需要额外系统依赖的包）。

5. **`3Dmol.js` 是一个 0 字节的空文件**，看起来是预留的分子可视化前端资源（3Dmol.js 是一个常用的网页端分子结构可视化库），但目前没有被任何脚本引用，也没有实际内容。如果暂时用不到，建议先删除，避免误导。

6. **没有自动化测试，也没有 CI 配置**，每一步是否跑对完全依赖人工检查终端输出和产出文件。

## 8. 背景知识速查（面向非化学背景读者）

- **EGFR / C797S**：EGFR 是驱动部分肺癌生长的关键蛋白，C797S 是其上一个氨基酸突变位点，会让第三代靶向药奥希替尼失去结合能力，是目前临床上重要的耐药机制之一。
- **IC50**：半数抑制浓度，数值越小代表化合物抑制靶点活性所需的浓度越低，即活性越强。
- **Lipinski's Rule of Five**：药物化学里判断一个小分子"是否具备成药潜力"的经验规则（分子量、脂水分配系数、氢键供受体数各有上限），不满足规则的分子通常口服吸收差，很少能成药。
- **Morgan 指纹 / Tanimoto 相似度**：把分子结构编码成一串二进制指纹（记录哪些局部原子环境出现过），再用 Tanimoto 系数比较两个指纹的重合程度，是化学信息学里最常用的"结构相似度"量化方法之一。
- **分子对接（Molecular Docking）/ AutoDock Vina**：把小分子"摆"进蛋白的结合口袋里，搜索最可能的结合姿势并给出一个能量打分（结合能，kcal/mol），分数越负通常代表预测的结合越稳定。这是一种计算预测，不等于真实实验结果，通常只用于化合物优先级排序、缩小后续湿实验验证的范围。
- **PDB / PDBQT**：PDB 是蛋白质/小分子三维结构的标准文本格式；PDBQT 是 AutoDock/Vina 专用的格式，在 PDB 基础上额外记录了原子类型、部分电荷和可旋转键信息，不能简单地把 PDB 文件改个后缀名当 PDBQT 用（对应本文档第 7 节的问题 1）。
