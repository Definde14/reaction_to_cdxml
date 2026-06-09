# reaction_to_cdxml

将反应路径（Reaction SMILES）转换为 **ChemDraw CDXML 文件**，可直接在 ChemDraw 中打开。

## 功能

- **单步/多步反应**：支持线性多步反应路径
- **反应条件标注**：箭头上方和下方添加条件文字
- **多组分反应**：自动处理多个反应物/产物，添加 "+" 分隔
- **电荷保留**：完整保留 `[N+]`、`[O-]`、`[Na+]` 等带电基团
- **文件批量输入**：从文件读取多步反应

## 依赖

```bash
pip install rdkit
```

## 使用方法

```bash
# 单步反应
python reaction_to_cdxml.py --rxn "c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]" -o nitration.cdxml

# 多步反应
python reaction_to_cdxml.py --rxns \
  "c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]" \
  "c1ccc(cc1)[N+](=O)[O-]>>c1ccc(cc1)N" \
  --labels "Step 1: Nitration" "Step 2: Reduction" \
  --conditions-above "HNO3/H2SO4" "H2/Pd" \
  --conditions-below "0 C to rt" "EtOH" \
  -o twostep.cdxml

# 从文件读取
python reaction_to_cdxml.py --file reactions.txt -o output.cdxml
```

## 输入格式

每步反应使用 Reaction SMILES 格式：

```
reactant1.reactant2>>product1.product2
```

例如硝化反应：
```
c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]
```

## 输出

生成标准的 `.cdxml` 文件，可以通过 ChemDraw 直接打开，支持编辑和导出为其他格式。