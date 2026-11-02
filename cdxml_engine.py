"""Reaction SMILES to editable CDXML; all geometry is measured in points."""
from dataclasses import dataclass
import math
import unicodedata
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import rdDepictor

DTD_URL = 'https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd'
BOND_LENGTH = 24.0
PAGE_WIDTH, PAGE_HEIGHT = 842.0, 595.0
MARGIN, GAP, ROW_GAP = 36.0, 24.0, 30.0


@dataclass
class ReactionStep:
    reactants_smiles: list[str]
    products_smiles: list[str]
    label: str = ''
    conditions_above: str = ''
    conditions_below: str = ''


def parse_reaction_smiles(rxn_smiles: str) -> tuple[list[str], list[str]]:
    parts = rxn_smiles.strip().split('>>')
    if len(parts) != 2 or not all(parts):
        raise ValueError('Expected nonempty reactants>>products')
    sides = [part.split('.') for part in parts]
    if any(not token.strip() for side in sides for token in side):
        raise ValueError('Empty molecule component')
    return tuple(sides)


def text_size(text: str, size: float = 10) -> tuple[float, float]:
    """Conservative font-independent envelope, including CJK and multiline text."""
    lines = text.split('\n')
    width = max(sum(1.1 if unicodedata.east_asian_width(c) in 'WF' else 0.85
                    for c in line) for line in lines) * size
    return max(size, width) + 8, len(lines) * size * 1.5 + 4


def generate_molecule_data(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError(f'Invalid SMILES: {smiles}')
    if any(a.HasQuery() or a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError('Query atoms and wildcard labels are not supported')
    if mol.GetStereoGroups():
        raise ValueError('Enhanced stereo groups are not supported')
    Chem.Kekulize(mol, clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    conf = mol.GetConformer()
    coords = [(conf.GetAtomPosition(i).x * BOND_LENGTH / 1.5,
               -conf.GetAtomPosition(i).y * BOND_LENGTH / 1.5)
              for i in range(mol.GetNumAtoms())]
    # Reserve labels on both sides, including isotope, charge and hydrogens.
    boxes = []
    for a, (x, y) in zip(mol.GetAtoms(), coords):
        label = a.GetSymbol() + 'H' + str(a.GetTotalNumHs())
        if a.GetIsotope():
            label += str(a.GetIsotope())
        if a.GetFormalCharge():
            label += str(a.GetFormalCharge())
        pad = text_size(label, 10)[0] / 2 if a.GetAtomicNum() != 6 or a.GetDegree() == 0 or a.GetIsotope() else 4
        boxes.append((x-pad, y-12, x+pad, y+12))
    left, top = min(b[0] for b in boxes), min(b[1] for b in boxes)
    right, bottom = max(b[2] for b in boxes), max(b[3] for b in boxes)
    return dict(mol=mol, coords=[(x-left,y-top) for x,y in coords],
                width=right-left, height=bottom-top)


class CDXMLBuilder:
    """XML serializer with globally unique IDs and explicit reaction semantics."""
    def __init__(self):
        self.next_id = 3
        self.root = ET.Element('CDXML', BondLength=str(BOND_LENGTH),
                               LabelFont='3', LabelSize='10', CaptionFont='3',
                               CaptionSize='10', LineWidth='0.6', BoldWidth='2')
        fonts = ET.SubElement(self.root, 'fonttable')
        ET.SubElement(fonts, 'font', id='3', charset='utf-8', name='Arial')

    def element(self, parent, tag: str, **attrs):
        self.next_id += 1
        return ET.SubElement(parent, tag, id=str(self.next_id),
                             **{k: str(v) for k,v in attrs.items()})

    def add_fragment(self, page, data: dict, x: float, y: float):
        frag = self.element(page, 'fragment', BoundingBox=box(x,y,x+data['width'],y+data['height']))
        ids = []
        for atom, (px,py) in zip(data['mol'].GetAtoms(), data['coords']):
            attrs = dict(p=f'{x+px:.3f} {y+py:.3f}', Element=atom.GetAtomicNum(),
                         NumHydrogens=atom.GetTotalNumHs(), Charge=atom.GetFormalCharge())
            if atom.GetIsotope():
                attrs['Isotope'] = atom.GetIsotope()
            radicals = atom.GetNumRadicalElectrons()
            if radicals:
                if radicals > 2:
                    raise ValueError('More than two radical electrons are not supported')
                attrs['Radical'] = 'Doublet' if radicals == 1 else 'Triplet'
            node = self.element(frag, 'n', **attrs)
            ids.append(node.get('id'))
        for bond in data['mol'].GetBonds():
            order = bond.GetBondTypeAsDouble()
            if order not in (1,2,3):
                raise ValueError(f'Unsupported bond type: {bond.GetBondType()}')
            attrs = dict(B=ids[bond.GetBeginAtomIdx()], E=ids[bond.GetEndAtomIdx()], Order=int(order))
            direction = bond.GetBondDir()
            if direction == Chem.BondDir.BEGINWEDGE:
                attrs['Display'] = 'WedgeBegin'
            elif direction == Chem.BondDir.BEGINDASH:
                attrs['Display'] = 'WedgedHashBegin'
            elif direction == Chem.BondDir.UNKNOWN:
                attrs['Display'] = 'Wavy'
            self.element(frag, 'b', **attrs)
        return frag.get('id')

    def add_text(self, page, text: str, x: float, y: float, size=10, bold=False):
        width, height = text_size(text,size)
        node = self.element(page, 't', p=f'{x:.3f} {y+size:.3f}',
                            BoundingBox=box(x,y,x+width,y+height), Justification='Left')
        ET.SubElement(node,'s', font='3', size=str(size), face='1' if bold else '0').text = text
        return node.get('id')

    def build(self) -> str:
        ET.indent(self.root)
        return '<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE CDXML SYSTEM "'+DTD_URL+'">\n'+ET.tostring(self.root,encoding='unicode')+'\n'


def box(*values: float) -> str:
    return ' '.join(f'{v:.3f}' for v in values)


def build_linear_reaction(steps: list[ReactionStep]) -> str:
    """Measure first, place left-to-right, then paginate without shrinking bonds.

    Wide rows use an expanded, integral A4 sheet grid; tall rows expand vertically.
    Whole reactions remain together. Labels occupy a dedicated band above the row.
    """
    if not steps:
        raise ValueError('At least one reaction is required')
    measured = []
    for step in steps:
        if not isinstance(step, ReactionStep):
            raise TypeError('Expected ReactionStep instances')
        if not step.reactants_smiles or not step.products_smiles:
            raise ValueError('Each reaction requires reactants and products')
        reactants = [generate_molecule_data(s) for s in step.reactants_smiles]
        products = [generate_molecule_data(s) for s in step.products_smiles]
        above, below = step.conditions_above, step.conditions_below
        arrow_w = max(80, *(text_size(t)[0] for t in (above,below)))
        mol_h = max(d['height'] for d in reactants+products)
        upper = max(mol_h/2, text_size(above)[1]+12 if above else 0)
        lower = max(mol_h/2, text_size(below)[1]+12 if below else 0)
        label_h = text_size(step.label,12)[1]+12 if step.label else 0
        widths = [sum(d['width'] for d in side) + max(0,len(side)-1)*40 for side in (reactants,products)]
        row_w = sum(widths)+2*GAP+arrow_w
        row_w = max(row_w,text_size(step.label,12)[0] if step.label else 0)
        measured.append((step,reactants,products,arrow_w,upper,lower,label_h,row_w))
    builder = CDXMLBuilder()
    page = None
    y = MARGIN
    width = math.ceil((max(row[-1] for row in measured)+2*MARGIN)/PAGE_WIDTH)*PAGE_WIDTH
    for step,reactants,products,arrow_w,upper,lower,label_h,row_w in measured:
        row_h = upper+lower+label_h
        if page is None or y+row_h > PAGE_HEIGHT-MARGIN:
            height = math.ceil((row_h+2*MARGIN)/PAGE_HEIGHT)*PAGE_HEIGHT
            page = builder.element(builder.root,'page', BoundingBox=box(0,0,width,height),
                                   WidthPages=int(width/PAGE_WIDTH),HeightPages=int(height/PAGE_HEIGHT))
            y = MARGIN
        center = y+label_h+upper
        if step.label:
            builder.add_text(page,step.label,MARGIN,y,12,True)
        x = MARGIN
        def side(items):
            nonlocal x
            result, plus = [], []
            for i,d in enumerate(items):
                if i:
                    plus.append(builder.add_text(page,'+',x+12,center-9,12))
                    x += 40
                result.append(builder.add_fragment(page,d,x,center-d['height']/2))
                x += d['width']
            return result,plus
        r_ids,r_plus = side(reactants)
        x += GAP
        arrow = builder.element(page,'arrow', Tail3D=f'{x:.3f} {center:.3f} 0',
                                Head3D=f'{x+arrow_w:.3f} {center:.3f} 0',
                                BoundingBox=box(x,center-3,x+arrow_w,center+3),
                                ArrowheadHead='Full',ArrowheadType='Solid',HeadSize='1000')
        above_ids,below_ids = [],[]
        if step.conditions_above:
            w,h = text_size(step.conditions_above)
            above_ids.append(builder.add_text(page,step.conditions_above,x+(arrow_w-w)/2,center-12-h))
        if step.conditions_below:
            w,h = text_size(step.conditions_below)
            below_ids.append(builder.add_text(page,step.conditions_below,x+(arrow_w-w)/2,center+12))
        x += arrow_w+GAP
        p_ids,p_plus = side(products)
        scheme = builder.element(page,'scheme')
        builder.element(scheme,'step',ReactionStepReactants=' '.join(r_ids),
                        ReactionStepProducts=' '.join(p_ids),ReactionStepArrows=arrow.get('id'),
                        ReactionStepPlusses=' '.join(r_plus+p_plus),
                        ReactionStepObjectsAboveArrow=' '.join(above_ids),
                        ReactionStepObjectsBelowArrow=' '.join(below_ids))
        y += row_h+ROW_GAP
    return builder.build()


