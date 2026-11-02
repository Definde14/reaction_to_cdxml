"""Real-paper regression: Tang et al., DOI 10.1002/anie.202110351, Fig. 3A."""
from collections import Counter
from itertools import combinations
from xml.etree import ElementTree as ET
import unittest
from rdkit import Chem
from cdxml_engine import build_linear_reaction
from paper_case_example import load_steps
import test_cdxml


class PaperCaseTests(unittest.TestCase):
    def test_paper_case(self) -> None:
        steps = load_steps()
        for step in steps:
            test_cdxml.ExportTests().document(
                ['.'.join(step.reactants_smiles) + '>>' + '.'.join(step.products_smiles)],
                label=step.label, conditions_above=step.conditions_above,
                conditions_below=step.conditions_below)
            sides = []
            for side in (step.reactants_smiles, step.products_smiles):
                counts = Counter()
                for smiles in side:
                    counts.update(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms())
                sides.append(counts)
            self.assertEqual(*sides)
        xml = build_linear_reaction(steps)
        expected = Counter(Chem.MolToSmiles(Chem.MolFromSmiles(s))
                           for step in steps for side in (step.reactants_smiles, step.products_smiles) for s in side)
        self.assertEqual(expected, Counter(Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(xml)))
        root = ET.fromstring(xml)
        ids = [e.get('id') for e in root.iter() if e.get('id')]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(len(root.findall('.//fragment')), 9)
        for page in root.findall('page'):
            _, _, width, height = map(float, page.get('BoundingBox').split())
            boxes = [tuple(map(float, e.get('BoundingBox').split())) for e in page if e.get('BoundingBox')]
            for left, top, right, bottom in boxes:
                self.assertTrue(0 <= left <= right <= width and 0 <= top <= bottom <= height)
            for a, b in combinations(boxes, 2):
                self.assertFalse(min(a[2], b[2])-max(a[0], b[0]) > .01 and
                                 min(a[3], b[3])-max(a[1], b[1]) > .01)


if __name__ == '__main__':
    unittest.main()

