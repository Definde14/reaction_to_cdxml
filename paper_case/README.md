# Tang 2021 论文案例

来源：[Tang et al., Angew. Chem. Int. Ed. 2021, 60, 23608–23613](https://doi.org/10.1002/anie.202110351)，Figure 3A 及第 23610 页正文。
手工转录前两步反应，不是 PDF 自动结构识别，也不包含完整的 HRP/TMB 显色步骤。

- 乳糖 + 水 → 葡萄糖 + 半乳糖，beta-Gal 催化。
- 葡萄糖 + 氧气 + 水 → 葡萄糖酸 + 过氧化氢：包含内酯水解的净反应。
- 水和副产物为配平补充；未添加未经原文核对的温度、时间、pH 或收率。
- 结构来自 PubChem 官方接口：乳糖 CID 440995、葡萄糖 5793、半乳糖 6036、葡萄糖酸 10690。
  API 模板：`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/property/IsomericSMILES,MolecularFormula/JSON`。
  还原端异头碳沿用数据库未指定状态，不代表单一异头体。

运行 `python paper_case_example.py -o tang2021_cascade.cdxml`，无需联网。
运行 `python -m unittest -v` 可复验几何、元素守恒和立体结构往返。

本地结果：9 个分子、78 个原子、74 条键、25 条楔线/虚楔线；画布为 1684 × 595 pt。
通过对象边界、包围盒不相交、唯一 ID、结构往返和元素守恒检查。
DTD 校验按项目 README 所述兼容处理官方 AS 声明。
该例需要两张横向 A4 的宽度，当前布局仍有较多留白。
RDKit 读取时出现 ReactionStepAgents fragment 查找警告，尽管全部分子可读取，
仍不能据此宣称反应语义在所有读取器中兼容。
ChemDraw COM 接口返回 0x80040112，原生显示与打印未验证。
