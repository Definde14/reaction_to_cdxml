"""Regression suite: geometry, official DTD, and chemical round trips."""
import os
import io
import re
import sys
import unittest
from pathlib import Path
from itertools import combinations
from xml.etree import ElementTree as ET
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rdkit import Chem
from cdxml_engine import ReactionStep, build_linear_reaction, parse_reaction_smiles


class ExportTests(unittest.TestCase):
    def document(self, reactions, **labels):
        steps = [ReactionStep(*parse_reaction_smiles(r), **labels) for r in reactions]
        xml = build_linear_reaction(steps)
        root = ET.fromstring(xml)
        ids = [e.get('id') for e in root.iter() if e.get('id')]
        self.assertEqual(len(ids),len(set(ids)))
        for page in root.findall('page'):
            _,_,width,height = map(float,page.get('BoundingBox').split())
            boxes=[]
            for child in page:
                if child.tag not in ('fragment','t','arrow'):
                    continue
                bounds=tuple(map(float,child.get('BoundingBox').split()))
                l,t,r,b=bounds
                self.assertTrue(0<=l<=r<=width and 0<=t<=b<=height, (child.tag,bounds))
                boxes.append((child.tag,bounds))
            for a,b in combinations(boxes,2):
                x,y=a[1],b[1]
                self.assertFalse(min(x[2],y[2])-max(x[0],y[0])>0.01 and
                                 min(x[3],y[3])-max(x[1],y[1])>0.01,(a,b))
            for arrow in page.findall('arrow'):
                self.assertGreater(float(arrow.get('Head3D').split()[0]),
                                   float(arrow.get('Tail3D').split()[0]))
        dtd_path=os.environ.get('CDXML_DTD')
        if dtd_path:
            from lxml import etree
            # Official DTD contains illegal XML NMTOKEN values + and ?.
            # Relax only AS enums so standard XML validators can parse it.
            source=Path(dtd_path).read_text()
            source=re.sub(r'AS \([^)]*\)', 'AS CDATA', source)
            dtd=etree.DTD(io.StringIO(source))
            self.assertTrue(dtd.validate(etree.fromstring(xml.encode())),str(dtd.error_log))
        return xml,root

    def test_chemistry_roundtrip(self):
        examples=['CCO','c1ccccc1[N+](=O)[O-]','[Na+]','[Cl-]',
                  '[13CH3]O','N[C@@H](C)C(=O)O','F/C=C/F','[CH3]','C']
        for smiles in examples:
            with self.subTest(smiles=smiles):
                xml,_=self.document([f'{smiles}>>{smiles}'])
                mols=Chem.MolsFromCDXML(xml)
                self.assertEqual(len(mols),2)
                expected=Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
                for mol in mols:
                    self.assertEqual(Chem.MolToSmiles(mol),expected)

    def test_many_steps_and_long_text(self):
        _,root=self.document(['CCO.[Na+].[Cl-]>>CC=O.O']*15,
                             label='Step label\n第二行',
                             conditions_above='Very long reaction conditions '*8+'<&>',
                             conditions_below='EtOH\n80 C\n24 hours')
        self.assertGreater(len(root.findall('page')),1)
        self.assertGreater(float(root.find('page').get('BoundingBox').split()[2]),842)

    def test_tall_molecule(self):
        self.document(['C1'+ 'C'*100+'1>>O'])

    def test_invalid_input(self):
        for reaction in ('','CCO>O','>>O','O>>','O..C>>CC'):
            with self.assertRaises(ValueError):
                parse_reaction_smiles(reaction)
        with self.assertRaises(ValueError):
            build_linear_reaction([])
        with self.assertRaises(ValueError):
            build_linear_reaction([ReactionStep(['invalid'],['O'])])


if __name__=='__main__':
    unittest.main()
