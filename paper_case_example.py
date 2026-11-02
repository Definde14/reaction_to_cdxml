"""Reproduce the two-step Tang 2021 paper case without network access."""
import argparse
import json
from pathlib import Path
from cdxml_engine import ReactionStep, build_linear_reaction

CASE_PATH = Path(__file__).resolve().parent / 'paper_case' / 'reactions.json'


def load_steps() -> list[ReactionStep]:
    return [ReactionStep(**record) for record in json.loads(CASE_PATH.read_text(encoding='utf-8'))]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output', type=Path, default=Path('tang2021_cascade.cdxml'))
    args = parser.parse_args()
    args.output.write_text(build_linear_reaction(load_steps()), encoding='utf-8')
    print(f'OK: {args.output}')


if __name__ == '__main__':
    main()
