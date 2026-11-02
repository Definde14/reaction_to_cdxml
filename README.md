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
## 布局重构

要求 Python 3.10+。安装：`pip install -r requirements.txt`。
原 CLI 参数和 `ReactionStep` / `build_linear_reaction` 入口保留。
内部实现位于 `cdxml_engine.py`，不再依赖 RDKit 可选的 ChemDraw 写出扩展。

- 坐标单位统一为 point，键长 24 pt；只使用二维描绘，不进行三维力场优化。
- 反应物在左、产物在右，箭头通过 `Tail3D` / `Head3D` 明确方向。
- 分子标签留白、条件文字和步骤标题参与尺寸计算；支持换行和 XML 特殊字符。
- 整步反应不拆分；超高分页，超宽扩展为整数 A4 横向纸张网格。
  扩展画布可能需要多张纸打印，并不等于缩放到单张 A4。
- 输出字体表、页面边界、唯一对象 ID 和 `scheme/step` 反应关联。
- 保留常见电荷、同位素、自由基、楔线立体结构及双键几何。
  目前仅支持 `reactants>>products`，不支持 agent 段、查询原子、通配符、增强立体组和配位键。
  原子映射编号不会作为输出标签保留。
- 文字尺寸采用保守字体无关估计；实际字体替换和复杂分子内部标签仍需 ChemDraw 人工复核。

## 验证与示例

```bash
python -m unittest -v
python gen_test_schemes.py --output-dir examples
```

测试覆盖对象边界、对象间相交、方向、唯一 ID、长条件、多组分、分页，
以及 RDKit 重新读取后的 canonical isomeric SMILES 对比。
设置 `CDXML_DTD` 为本地官方 DTD 路径可额外运行格式验证（需 `pip install lxml`）。
官方 DTD 的 `AS` 枚举含 XML NMTOKEN 不允许的 `+` 和 `?`，测试仅将这一属性
声明放宽为 CDATA 后加载，其余声明保持原样。CI 自动下载 DTD 执行检查。
这些验证不替代 ChemDraw 原生打开和打印检查。

设计参考：

- [Revvity 官方 CDXML DTD](https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd)
- [RDKit CDXML 文档](https://rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html)
- [pycdxml](https://github.com/kienerj/pycdxml)：借鉴其统一样式和布局思路，未复制代码。
