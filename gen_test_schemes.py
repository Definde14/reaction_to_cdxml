"""Generate portable regression examples using the production exporter."""
import argparse
from pathlib import Path
from reaction_to_cdxml import ReactionStep, build_linear_reaction, parse_reaction_smiles


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',default='examples')
    args=parser.parse_args()
    out=Path(args.output_dir)
    out.mkdir(parents=True,exist_ok=True)
    examples={
        'nitration':[ReactionStep(*parse_reaction_smiles('c1ccccc1>>c1ccc(cc1)[N+](=O)[O-]'),
                                  'Nitration','HNO3 / H2SO4','0 C to rt')],
        'multicomponent':[ReactionStep(*parse_reaction_smiles('CCO.[Na+].[Cl-]>>CC=O.O'),
                                       'Oxidation','Catalyst, solvent\n24 hours','80 C')],
        'pagination':[ReactionStep(['CCO'],['CC=O'],f'Step {i+1}') for i in range(12)],
        'stereo':[ReactionStep(['N[C@@H](C)C(=O)O'],['N[C@@H](C)C(=O)OC'])],
    }
    for name,steps in examples.items():
        (out/f'{name}.cdxml').write_text(build_linear_reaction(steps),encoding='utf-8')


if __name__=='__main__':
    main()
