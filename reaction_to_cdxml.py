"""Command line interface and backward-compatible public imports."""
import argparse
import sys
from cdxml_engine import (ReactionStep, CDXMLBuilder, generate_molecule_data,
                          parse_reaction_smiles, build_linear_reaction)

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

    try:
        steps = []
        for i, rxn_smi in enumerate(rxns):
            reactants, products = parse_reaction_smiles(rxn_smi)
            label = args.labels[i] if args.labels and i < len(args.labels) else ""
            above = args.conditions_above[i] if i < len(args.conditions_above) else ""
            below = args.conditions_below[i] if i < len(args.conditions_below) else ""
            steps.append(ReactionStep(reactants, products, label, above, below))
        cdxml = build_linear_reaction(steps)
    except (ValueError, TypeError) as exc:
        parser.error(str(exc))

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
