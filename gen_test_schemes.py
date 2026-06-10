#!/usr/bin/env python3
"""Generate CDXML test files with different coordinate schemes."""
from rdkit import Chem
from rdkit.Chem import AllChem
import re, os

OUT = "d:/Coding/reaction_to_cdxml"

def parse_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    raw = Chem.MolToCDXMLBlock(mol)
    atoms = []
    for m in re.finditer(r'<n\s+id="(\d+)"\s+p="(-?[\d\.]+)\s+(-?[\d\.]+)"(.*?)/>', raw, re.DOTALL):
        atoms.append((int(m.group(1)), float(m.group(2)), float(m.group(3)), m.group(4)))
    bonds = []
    for m in re.finditer(r'<b\s+id="(\d+)"\s+B="(\d+)"\s+E="(\d+)"(.*?)/>', raw, re.DOTALL):
        bonds.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4)))
    xs = [a[1] for a in atoms]
    ys = [a[2] for a in atoms]
    return {
        'atoms': atoms, 'bonds': bonds,
        'minx': min(xs), 'maxx': max(xs), 'miny': min(ys), 'maxy': max(ys),
        'cx': (min(xs)+max(xs))/2, 'cy': (min(ys)+max(ys))/2,
        'w': max(xs)-min(xs), 'h': max(ys)-min(ys),
    }

def shift(mol, dx, dy):
    natoms = [(a[0], a[1]+dx, a[2]+dy, a[3]) for a in mol['atoms']]
    return natoms, mol['bonds']

def build(fragments, arrows=None, texts=None):
    nid = 100
    seg = []
    seg += ['<?xml version="1.0" encoding="UTF-8" ?>',
            '<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >',
            '<CDXML', ' id="1"', ' BondLength=""', '>',
            '<page', ' id="2">']
    for atoms, bonds in fragments:
        seg += ['<fragment', f' id="{nid}"', '>']; nid += 1
        id_map = {}
        for a in atoms: id_map[a[0]] = nid; nid += 1
        for a in atoms:
            seg += ['<n', f' id="{id_map[a[0]]}"', f' p="{a[1]:.2f} {a[2]:.2f}"']
            cleaned = ' '.join(a[3].split())
            if cleaned: seg += [f' {cleaned}']
            seg += ['/>']
        for b in bonds:
            seg += ['<b', f' id="{nid}"', f' B="{id_map[b[1]]}"', f' E="{id_map[b[2]]}"']
            cleaned = ' '.join(b[3].split())
            if cleaned: seg += [f' {cleaned}']
            seg += ['/>']; nid += 1
        seg += ['</fragment>']
    if arrows:
        for a in arrows:
            seg += ['<graphic', f' id="{nid}"', ' GraphicType="Line"',
                    f' BoundingBox="{a[0]:.2f} {a[1]:.2f} {a[2]:.2f} {a[3]:.2f}"',
                    ' ArrowType="FullHead"', '/>']; nid += 1
    if texts:
        for t in texts:
            seg += ['<t', f' id="{nid}"', f' p="{t[0]:.2f} {t[1]:.2f}"',
                    ' Justification="1">', '<s', ' font="2"', ' size="10"', ' face="0">',
                    t[2], '</s>', '</t>']; nid += 1
    seg += ['</page>', '</CDXML>']
    return '\r\n'.join(seg) + '\r\n'

# ---- 准备分子 ----
b = parse_mol('c1ccccc1')
n = parse_mol('c1ccc(cc1)[N+](=O)[O-]')
a = parse_mol('c1ccc(cc1)N')
gap = 50

# ---- A: RDKit原始风格（x有负值） ----
fs = [shift(b, 0, 0), shift(n, b['w']+gap, 0)]
ar = [(b['w']+gap//2, 0, b['w']+gap//2+60, 0)]
tx = [(-20, 20, "A: RDKit原始坐标(x有正有负)")]
with open(os.path.join(OUT, "test_A.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
# 检查方向
xs = [a[1] for a in fs[0][0] + fs[1][0]]
ys = [a[2] for a in fs[0][0] + fs[1][0]]
print(f"A: x=[{min(xs):.0f},{max(xs):.0f}] y=[{min(ys):.0f},{max(ys):.0f}] dir={'L->R' if max(xs)>0 else '?'}")

# ---- B: 全正x，单行 ----
fs = [shift(b, 0, 0), shift(n, b['w']+gap, 0)]
# 平移使最左=0
x0 = min(min(a[1] for a in fs[0][0]), min(a[1] for a in fs[1][0]))
fs = [([(a[0], a[1]-x0, a[2], a[3]) for a in fs[0][0]], fs[0][1]),
      ([(a[0], a[1]-x0, a[2], a[3]) for a in fs[1][0]], fs[1][1])]
xs = [a[1] for a in fs[0][0] + fs[1][0]]
ar_x = (max(a[1] for a in fs[0][0]) + gap//2)
ar = [(ar_x, 0, ar_x+60, 0)]
tx = [(5, 20, f"B: 全正x x=[{min(xs):.0f},{max(xs):.0f}]")]
with open(os.path.join(OUT, "test_B.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
print(f"B: x=[{min(xs):.0f},{max(xs):.0f}]")

# ---- C: 两行，Step1在上(y=80), Step2在下(y=-80) ----
fs = [shift(b, 0, 80), shift(n, b['w']+gap, 80),
      shift(n, 0, -80), shift(a, n['w']+gap, -80)]
# 统一平移使最左=0
all_xs = [a[1] for fr in fs for a in fr[0]]
x0 = min(all_xs)
fs = [([(a[0], a[1]-x0, a[2], a[3]) for a in fr[0]], fr[1]) for fr in fs]
ar = [(b['w']+gap//2, 80, b['w']+gap//2+60, 80),
      (n['w']+gap//2, -80, n['w']+gap//2+60, -80)]
tx = [(10, 105, "C: Step1 y=+80 (top), Step2 y=-80 (bottom)")]
with open(os.path.join(OUT, "test_C.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
ys = [a[2] for fr in fs for a in fr[0]]
print(f"C: x=[{min([a[1] for fr in fs for a in fr[0]]):.0f},{max([a[1] for fr in fs for a in fr[0]]):.0f}] y=[{min(ys):.0f},{max(ys):.0f}]")

# ---- D: Step1在负y, Step2在正y(y翻转假设) ----
fs = [shift(b, 0, -80), shift(n, b['w']+gap, -80),
      shift(n, 0, 80), shift(a, n['w']+gap, 80)]
x0 = min([a[1] for fr in fs for a in fr[0]])
fs = [([(a[0], a[1]-x0, a[2], a[3]) for a in fr[0]], fr[1]) for fr in fs]
ar = [(b['w']+gap//2, -80, b['w']+gap//2+60, -80),
      (n['w']+gap//2, 80, n['w']+gap//2+60, 80)]
tx = [(10, -105, "D: Step1 y=-80, Step2 y=+80 (y-flip)")]
with open(os.path.join(OUT, "test_D.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
ys = [a[2] for fr in fs for a in fr[0]]
print(f"D: y=[{min(ys):.0f},{max(ys):.0f}]")

# ---- E: 全正y，Step1在y=250, Step2在y=80 ----
fs = [shift(b, 0, 250), shift(n, b['w']+gap, 250),
      shift(n, 0, 80), shift(a, n['w']+gap, 80)]
x0 = min([a[1] for fr in fs for a in fr[0]])
fs = [([(a[0], a[1]-x0, a[2], a[3]) for a in fr[0]], fr[1]) for fr in fs]
ar = [(b['w']+gap//2, 250, b['w']+gap//2+60, 250),
      (n['w']+gap//2, 80, n['w']+gap//2+60, 80)]
tx = [(10, 275, "E: All y>0, Step1 y=250, Step2 y=80")]
with open(os.path.join(OUT, "test_E.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
ys = [a[2] for fr in fs for a in fr[0]]
print(f"E: y=[{min(ys):.0f},{max(ys):.0f}]")

# ---- F: 当前工具风格，Step1在y=0, Step2在y=-130 ----
fs = [shift(b, 0, 0), shift(n, b['w']+gap, 0),
      shift(n, 0, -130), shift(a, n['w']+gap, -130)]
x0 = min([a[1] for fr in fs for a in fr[0]])
fs = [([(a[0], a[1]-x0, a[2], a[3]) for a in fr[0]], fr[1]) for fr in fs]
ar = [(b['w']+gap//2, 0, b['w']+gap//2+60, 0),
      (n['w']+gap//2, -130, n['w']+gap//2+60, -130)]
tx = [(10, 25, "F: Step1 y=0, Step2 y=-130 (current tool style)")]
with open(os.path.join(OUT, "test_F.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
ys = [a[2] for fr in fs for a in fr[0]]
print(f"F: y=[{min(ys):.0f},{max(ys):.0f}]")

# ---- G: 全正y，Step1 在上(大y), Step2 在下(小y), 第一行用 y=400 ----
fs = [shift(b, 0, 400), shift(n, b['w']+gap, 400),
      shift(n, 0, 180), shift(a, n['w']+gap, 180)]
x0 = min([a[1] for fr in fs for a in fr[0]])
fs = [([(a[0], a[1]-x0+b['w']//3, a[2], a[3]) for a in fr[0]], fr[1]) for fr in fs]
ar = [(b['w']+gap//2, 400, b['w']+gap//2+60, 400),
      (n['w']+gap//2, 180, n['w']+gap//2+60, 180)]
tx = [(10, 425, "G: All y>0, Step1 y=400, Step2 y=180")]
with open(os.path.join(OUT, "test_G.cdxml"), 'wb') as f:
    f.write(build(fs, ar, tx).encode('utf-8'))
ys = [a[2] for fr in fs for a in fr[0]]
print(f"G: y=[{min(ys):.0f},{max(ys):.0f}]")

print()
print("请在 ChemDraw 中依次打开 test_A ~ test_G")
print("告诉我哪个方案的方向是 从左到右、从上到下 !")