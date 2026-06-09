#!/usr/bin/env python3
"""
reaction_to_cdxml.py
====================
将反应路径（Reaction SMILES 或序列）转换为 ChemDraw CDXML 文件，
可直接在 ChemDraw 中打开。

用法示例:
  # 单步反应: 苯 -> 硝基苯
  python reaction_to_cdxml.py --rxn "c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]" -o benzene_nitration.cdxml

  # 多步反应 (用空格分隔多个 SMILES)
  python reaction_to_cdxml.py --rxns \
    "c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]" \
    "c1ccc(cc1)[N+](=O)[O-]>>c1ccc(cc1)N" \
    -o twostep.cdxml

  # 指定步骤标签
  python reaction_to_cdxml.py --rxns "RXN_SMILES1" "RXN_SMILES2" \
    --labels "Step 1: Nitration" "Step 2: Reduction" \
    -o labeled.cdxml

  # 从文件读取反应 SMILES（每行一个）
  python reaction_to_cdxml.py --file reactions.txt -o output.cdxml

注意事项:
  - 反应 SMILES 格式:  reactant1.reactant2>>product1.product2
  - 多步反应: 用空格分隔多个 Reaction SMILES
  - 每个分子必须有有效 SMILES, RDKit 会自动生成 2D 坐标
  - 默认输出到当前目录的 reaction.cdxml

依赖:
  pip install rdkit
"""

import re
import sys
import argparse
from xml.sax.saxutils import escape as xml_escape
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


# ============================================================
# CDXML 常量
# ============================================================
DEFAULT_BOND_LENGTH = 1.5  # ChemDraw 默认键长 (cm)
SCALE_FACTOR = 28.8        # RDKit internal -> CDXML 坐标转换因子
                            # RDKit 中键长 ≈ 1.5 → CDXML 键长 ≈ 43.2
                            # CDXML 中两个相邻原子约 43.2 单位
# 布局参数
MOL_GAP = 80.0             # 分子之间的间距（≈2.8 个键长）
ARROW_LENGTH = 86.0        # 箭头长度（≈2 cm）
ARROW_GAP = 36.0           # 箭头与分子的间距（≈1.2 cm）
VERTICAL_OFFSET = 0.0      # 垂直居中偏移

# 页面设置（A4 横向: 297mm × 210mm）
# CDXML 单位: 1 inch = 72 units, 1 mm ≈ 2.835 units
PAGE_WIDTH = 842.0    # 297mm
PAGE_HEIGHT = 595.0   # 210mm

# 步骤标签字体设置
LABEL_FONT = 2    # Times
LABEL_SIZE = 12   # point


# ============================================================
# 核心: 解析 RDKit 生成的 CDXML, 提取原子数据
# ============================================================

def parse_cdxml_atoms(cdxml_text: str):
    """解析 RDKit 生成的 CDXML, 提取所有原子 (节点) 信息。

    返回:
        atoms: list of dict, 每个 dict 有 keys:
            id: int, 节点 ID
            x: float, y: float, 原始坐标
            element: int, 原子序数 (0=未指定, 通常为 C)
            raw_xml: str, 原始 <n .../> 字符串
        bonds: list of dict, 每个 dict 有 keys:
            id: int
            B: int, Begin 节点 ID
            E: int, End 节点 ID
            Order: float or None
            raw_xml: str
    """
    atoms = []
    bonds = []

    # 用正则解析节点 - 匹配所有可能的属性，确保不丢失任何节点
    # CDXML 中的 <n> 元素可以有 id, p, Element, NumHydrogens, Charge, Radical, Isotope 等
    # 我们需要能够捕获所有这些属性的行
    node_pat = re.compile(r'<n\s+id="(\d+)"\s+p="(-?[\d\.]+)\s+(-?[\d\.]+)"(.*?)/>', re.DOTALL)
    for match in node_pat.finditer(cdxml_text):
        node_id = int(match.group(1))
        x = float(match.group(2))
        y = float(match.group(3))
        rest = match.group(4)
        # 从剩余字符串中提取 Element 属性（如果有）
        elem_match = re.search(r'Element="(\d+)"', rest)
        element = int(elem_match.group(1)) if elem_match else 6
        atoms.append({
            'id': node_id,
            'x': x,
            'y': y,
            'element': element,
            'attrs': rest.strip(),  # 保存其他属性用于重建
        })

    # 用正则解析键 - 匹配所有可能的属性
    bond_pat = re.compile(r'<b\s+id="(\d+)"\s+B="(\d+)"\s+E="(\d+)"(.*?)/>', re.DOTALL)
    for match in bond_pat.finditer(cdxml_text):
        bond_id = int(match.group(1))
        B = int(match.group(2))
        E = int(match.group(3))
        rest = match.group(4)
        order_match = re.search(r'Order="([\d\.]+)"', rest)
        Order = float(order_match.group(1)) if order_match else None
        bonds.append({
            'id': bond_id,
            'B': B,
            'E': E,
            'Order': Order,
            'attrs': rest.strip(),
        })

    return atoms, bonds


def get_molecule_center(atoms):
    """获取分子几何中心"""
    if not atoms:
        return 0, 0
    xs = [a['x'] for a in atoms]
    ys = [a['y'] for a in atoms]
    return (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2


# ============================================================
# 生成单分子 CDXML
# ============================================================

def generate_molecule_data(smiles: str) -> tuple:
    """生成分子数据（坐标、CDXML），不放置到具体位置。

    返回:
        (atoms, bounds, width, height, center_x, center_y)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无效 SMILES: {smiles}")

    AllChem.Compute2DCoords(mol)
    if mol.GetNumAtoms() > 1:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except:
            pass

    cdxml = Chem.MolToCDXMLBlock(mol)
    atoms, bonds = parse_cdxml_atoms(cdxml)

    if not atoms:
        raise ValueError(f"分子 {smiles} 没有可解析的原子")

    cx, cy = get_molecule_center(atoms)
    xs = [a['x'] for a in atoms]
    ys = [a['y'] for a in atoms]
    width = max(xs) - min(xs)
    height = max(ys) - min(ys)

    return {
        'atoms': atoms,
        'bonds': bonds,
        'center_x': cx,
        'center_y': cy,
        'width': width,
        'height': height,
    }


# ============================================================
# 核心: 多步反应 CDXML 生成器
# ============================================================

class ReactionStep:
    """一个反应步骤"""
    def __init__(self, reactants_smiles: list, products_smiles: list,
                 label: str = "", conditions_above: str = "",
                 conditions_below: str = ""):
        self.reactants_smiles = reactants_smiles
        self.products_smiles = products_smiles
        self.label = label
        self.conditions_above = conditions_above
        self.conditions_below = conditions_below


class CDXMLBuilder:
    """CDXML 多步反应布局构建器"""

    def __init__(self):
        self.next_id = 1
        self.fragments = []
        self.graphics = []
        self.texts = []
        self.arrows = []

    def _get_id(self):
        _id = self.next_id
        self.next_id += 1
        return _id

    def add_fragment(self, mol_data: dict, offset_x: float = 0.0,
                  offset_y: float = 0.0) -> tuple:
        """基于预生成的分子数据添加 fragment。

        Args:
            mol_data: generate_molecule_data() 返回的 dict
            offset_x, offset_y: 放置位置偏移（从分子几何中心计算）

        返回: (fragment_id, bounds)
        """
        atoms = []
        for a in mol_data['atoms']:
            atoms.append({
                'id': a['id'],
                'x': a['x'],
                'y': a['y'],
                'element': a['element'],
                'attrs': a.get('attrs', ''),
            })
        bonds = [dict(b) for b in mol_data['bonds']]

        # 计算偏移：放置位置的几何中心
        cx, cy = mol_data['center_x'], mol_data['center_y']
        dx = offset_x - cx
        dy = offset_y - cy

        # 重新映射所有 ID 为全局唯一
        id_map = {}
        for atom in atoms:
            old_id = atom['id']
            new_id = self._get_id()
            id_map[old_id] = new_id
            atom['id'] = new_id
            atom['x'] += dx
            atom['y'] += dy
            atom['attrs'] = ' '.join(atom.get('attrs', '').split())

        for bond in bonds:
            bond['id'] = self._get_id()
            bond['B'] = id_map[bond['B']]
            bond['E'] = id_map[bond['E']]

        fragment_id = self._get_id()
        self.fragments.append({
            'id': fragment_id,
            'atoms': atoms,
            'bonds': bonds,
        })

        return fragment_id

    def add_arrow(self, x1: float, y1: float, x2: float, y2: float,
                  arrow_type: str = "FullHead"):
        """添加反应箭头"""
        arrow_id = self._get_id()
        self.arrows.append({
            'id': arrow_id,
            'x1': x1, 'y1': y1,
            'x2': x2, 'y2': y2,
            'arrow_type': arrow_type,
        })
        return arrow_id

    def add_text(self, text: str, x: float, y: float,
                 font_size: int = LABEL_SIZE, font_face: int = LABEL_FONT,
                 justification: int = 0, bold: bool = False):
        """添加文本标签"""
        text_id = self._get_id()
        self.texts.append({
            'id': text_id,
            'text': text,
            'x': x, 'y': y,
            'font_size': font_size,
            'font_face': font_face,
            'justification': justification,
            'bold': bold,
        })
        return text_id

    def add_plus(self, x: float, y: float):
        """添加加号 (组分分隔符)"""
        plus_id = self._get_id()
        self.texts.append({
            'id': plus_id,
            'text': '+',
            'x': x, 'y': y,
            'font_size': 18,
            'font_face': 2,
            'justification': 1,  # center
            'bold': False,
        })
        return plus_id

    def build(self) -> str:
        """生成最终 CDXML 字符串（与 RDKit 格式一致）"""
        L = []  # each entry = one line, no trailing CRLF

        L.append('<?xml version="1.0" encoding="UTF-8" ?>')
        L.append('<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >')

        L.append('<CDXML')
        L.append(f' id="{self._get_id()}"')
        L.append(' BondLength=""')
        L.append('>')

        L.append('<page')
        L.append(f' id="{self._get_id()}">')

        for frag in self.fragments:
            L.append('<fragment')
            L.append(f' id="{frag["id"]}">')

            for atom in frag['atoms']:
                x_str = f"{atom['x']:.2f}"
                y_str = f"{atom['y']:.2f}"
                extra = atom.get('attrs', '')
                L.append('<n')
                L.append(f' id="{atom["id"]}"')
                L.append(f' p="{x_str} {y_str}"')
                if extra:
                    L.append(f' {extra}')
                elif atom['element'] != 6:
                    L.append(f' Element="{atom["element"]}"')
                L.append('/>')

            for bond in frag['bonds']:
                extra = bond.get('attrs', '')
                L.append('<b')
                L.append(f' id="{bond["id"]}"')
                L.append(f' B="{bond["B"]}"')
                L.append(f' E="{bond["E"]}"')
                if extra:
                    L.append(f' {extra}')
                L.append('/>')

            L.append('</fragment>')

        for arrow in self.arrows:
            L.append('<graphic')
            L.append(f' id="{arrow["id"]}"')
            L.append(' GraphicType="Line"')
            L.append(f' BoundingBox="{arrow["x1"]:.2f} {arrow["y1"]:.2f} {arrow["x2"]:.2f} {arrow["y2"]:.2f}"')
            L.append(f' ArrowType="{arrow["arrow_type"]}"')
            L.append('/>')

        for text_item in self.texts:
            L.append('<t')
            L.append(f' id="{text_item["id"]}"')
            L.append(f' p="{text_item["x"]:.2f} {text_item["y"]:.2f}"')
            L.append(f' Justification="{text_item["justification"]}">')
            L.append('<s')
            L.append(f' font="{text_item["font_face"]}"')
            L.append(f' size="{text_item["font_size"]}"')
            L.append(f' face="{1 if text_item["bold"] else 0}">')
            L.append(f'{xml_escape(text_item["text"])}')
            L.append('</s>')
            L.append('</t>')

        L.append('</page>')
        L.append('</CDXML>')

        return '\r\n'.join(L) + '\r\n'


def build_linear_reaction(steps: list) -> str:
    """构建线性多步反应 CDXML。

    布局（每行居中）:
      Reactant1 + Reactant2 → Product1 + Product2   (Step 1)
                ↓
      Reactant3 → Product3                          (Step 2)

    Args:
        steps: list of ReactionStep 对象

    Returns:
        CDXML 字符串
    """
    builder = CDXMLBuilder()

    # 标准化输入
    processed_steps = []
    for i, step in enumerate(steps):
        if isinstance(step, ReactionStep):
            processed_steps.append(step)
        else:
            raise TypeError(f"steps[{i}] 格式无效")

    # Phase 1: 预计算所有分子的数据 (几何中心, 宽, 高)
    step_mol_data = []
    for step in processed_steps:
        r_data = [generate_molecule_data(smi) for smi in step.reactants_smiles]
        p_data = [generate_molecule_data(smi) for smi in step.products_smiles]
        step_mol_data.append((r_data, p_data))

    # Phase 2: 布局计算
    # 每行从左到右:
    #   R1 | + | R2 | → | P1 | + | P2
    # 间距常量: MOL_GAP 是分子间的加号间距
    #           分子间间距 = gap_rp (分子-加号-分子)

    current_y = 0.0  # 行中心 y 坐标（ChemDraw 坐标系：从下到上为正）

    for step_idx, (step, (r_data, p_data)) in enumerate(zip(processed_steps, step_mol_data)):
        # 计算该行最大高度（用于行间距）
        max_h = 0.0
        for d in r_data + p_data:
            if d['height'] > max_h:
                max_h = d['height']

        y_center = current_y

        # ========= 布局计算 =========
        # 思路：先确定所有元素（分子、加号、箭头）的左中右位置
        # 每个分子用它的 center_x 定位；箭头用 x1,x2 定位

        # ---- 反应物区 ----
        rx = 0.0  # 当前 x 游标
        r_centers = []   # 每个反应物几何中心应放置的 x 坐标
        r_rights = []     # 每个反应物右边界
        for d in r_data:
            half_w = d['width'] / 2
            cx = rx + half_w  # 分子的几何中心
            r_centers.append(cx)
            r_rights.append(cx + half_w)
            rx = cx + half_w + MOL_GAP * 0.4  # 加号间距

        # ---- 反应物之间的加号 ----
        r_plus_positions = []
        if len(r_data) > 1:
            for i in range(len(r_data) - 1):
                # 加号放在两个分子之间的中点
                plus_x = (r_rights[i] + r_centers[i+1] - r_data[i+1]['width'] / 2) / 2
                r_plus_positions.append(plus_x)

        # ---- 箭头 ----
        if r_data:
            arrow_start = r_rights[-1] + ARROW_GAP
        else:
            arrow_start = rx + ARROW_GAP
        arrow_end = arrow_start + ARROW_LENGTH
        arrow_cx = (arrow_start + arrow_end) / 2

        # ---- 产物区 ----
        px = arrow_end + ARROW_GAP
        p_centers = []   # 每个产物几何中心应放置的 x 坐标
        p_rights = []     # 每个产物右边界
        for d in p_data:
            half_w = d['width'] / 2
            cx = px + half_w
            p_centers.append(cx)
            p_rights.append(cx + half_w)
            px = cx + half_w + MOL_GAP * 0.4

        # ---- 产物之间的加号 ----
        p_plus_positions = []
        if len(p_data) > 1:
            for i in range(len(p_data) - 1):
                plus_x = (p_rights[i] + p_centers[i+1] - p_data[i+1]['width'] / 2) / 2
                p_plus_positions.append(plus_x)

        # ========= 输出 =========
        # 放置反应物
        for d, cx in zip(r_data, r_centers):
            builder.add_fragment(d, cx, y_center)

        # 放置反应物加号
        for plus_x in r_plus_positions:
            builder.add_plus(plus_x, y_center)

        # 箭头
        builder.add_arrow(arrow_start, y_center, arrow_end, y_center)

        # 放置产物
        for d, cx in zip(p_data, p_centers):
            builder.add_fragment(d, cx, y_center)

        # 放置产物加号
        for plus_x in p_plus_positions:
            builder.add_plus(plus_x, y_center)

        # 条件文字
        step_obj = processed_steps[step_idx]
        if step_obj.conditions_above:
            builder.add_text(step_obj.conditions_above, arrow_cx, y_center + 14,
                             font_size=10, justification=1)
        if step_obj.conditions_below:
            builder.add_text(step_obj.conditions_below, arrow_cx, y_center - 14,
                             font_size=10, justification=1)
        if step_obj.label:
            builder.add_text(step_obj.label, arrow_cx,
                             y_center + max_h / 2 + 20,
                             font_size=12, font_face=2, justification=1, bold=True)

        # 行间距（向负y方向移动，第一行在顶部大y，后续行在下方小y）
        # 但在 ChemDraw 坐标系中 y > 0，所以需确保所有 y 正值
        # 最终 builder.build() 前会整体平移使 min_y >= 页面底部
        current_y -= (max_h + MOL_GAP * 1.2)

    return builder.build()


# ============================================================
# 解析反应 SMILES
# ============================================================

def parse_reaction_smiles(rxn_smiles: str):
    """解析反应 SMILES, 返回 (reactants, products) 列表"""
    # 先用 RDKit 的 reaction SMILES 解析
    parts = rxn_smiles.split('>>')
    if len(parts) != 2:
        raise ValueError(f"无效反应 SMILES (需要包含 >>): {rxn_smiles}")

    reactant_part, product_part = parts
    reactants = reactant_part.split('.') if reactant_part else []
    products = product_part.split('.') if product_part else []

    if not reactants:
        raise ValueError(f"反应 SMILES 缺少反应物: {rxn_smiles}")
    if not products:
        raise ValueError(f"反应 SMILES 缺少产物: {rxn_smiles}")

    return reactants, products


def read_reactions_from_file(filepath: str):
    """从文件读取反应 (每行一个 SMILES, 空行和 # 开头为注释)"""
    steps = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            steps.append(line)
    return steps


# ============================================================
# CLI 入口
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='将反应路径 (Reaction SMILES) 转换为 ChemDraw CDXML 文件',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s --rxn "c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]" -o nitration.cdxml
  %(prog)s --rxns "A>>B" "B>>C" -o twostep.cdxml
  %(prog)s --rxns "A>>B" --labels "Step 1" --conditions-above "HNO3/H2SO4" -o labeled.cdxml
  %(prog)s --file reactions.txt -o output.cdxml
        """
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--rxn', type=str,
                       help='单步反应 SMILES (格式: reactants>>products)')
    group.add_argument('--rxns', type=str, nargs='+',
                       help='多步反应 SMILES (多个 SMILES 用空格分隔)')
    group.add_argument('--file', type=str,
                       help='从文件读取反应 (每行一个 SMILES, # 开头为注释)')

    parser.add_argument('-o', '--output', type=str, default='reaction.cdxml',
                        help='输出 CDXML 文件路径 (默认: reaction.cdxml)')
    parser.add_argument('--labels', type=str, nargs='*',
                        help='步骤标签 (如 "Step 1" "Step 2")')
    parser.add_argument('--conditions-above', type=str, nargs='*', default=[],
                        help='箭头上方条件文字')
    parser.add_argument('--conditions-below', type=str, nargs='*', default=[],
                        help='箭头下方条件文字')

    args = parser.parse_args()

    # 收集反应步骤
    if args.rxn:
        rxns = [args.rxn]
    elif args.rxns:
        rxns = args.rxns
    elif args.file:
        rxns = read_reactions_from_file(args.file)
    else:
        parser.error("需要提供 --rxn, --rxns 或 --file")

    steps = []
    for i, rxn_smi in enumerate(rxns):
        reactants, products = parse_reaction_smiles(rxn_smi)
        label = args.labels[i] if args.labels and i < len(args.labels) else ""
        above = args.conditions_above[i] if args.conditions_above and i < len(args.conditions_above) else ""
        below = args.conditions_below[i] if args.conditions_below and i < len(args.conditions_below) else ""
        steps.append(ReactionStep(reactants, products, label, above, below))

    # 生成 CDXML
    try:
        cdxml = build_linear_reaction(steps)
    except Exception as e:
        sys.stderr.reconfigure(encoding='utf-8')
        print(f"CDXML ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # 写入文件（用二进制模式防止 Windows 自动转换 \n → \r\n）
    with open(args.output, 'wb') as f:
        f.write(cdxml.encode('utf-8'))

    sys.stdout.reconfigure(encoding='utf-8')  # for Windows GBK
    print(f"  OK: {args.output}")
    print(f"  反应步数: {len(steps)}")
    for i, step in enumerate(steps):
        r_str = ' + '.join(step.reactants_smiles)
        p_str = ' + '.join(step.products_smiles)
        print(f"  {i+1}. {r_str} -> {p_str}")


if __name__ == '__main__':
    main()