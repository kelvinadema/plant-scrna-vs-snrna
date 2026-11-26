#!/usr/bin/env python3
"""Run the current Chapter 3 analyses from frozen processed inputs."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section', nargs='+', help='Section numbers, e.g. 04 08; default: all')
    parser.add_argument('--figures', action='store_true', help='Also render figures and supplementary tables with R')
    parser.add_argument('--figures-only', action='store_true', help='Render previously generated analysis tables')
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[key] = '1'
    os.environ['MPLBACKEND'] = 'Agg'
    sources = sorted((ROOT/'scripts/analysis').glob('[0-9][0-9]_*.py'))
    if args.section:
        selected = {number.zfill(2) for number in args.section}
        unknown = selected - {source.name[:2] for source in sources}
        if unknown:
            parser.error(f'Unknown sections: {sorted(unknown)}')
        sources = [source for source in sources if source.name[:2] in selected]
    records = []
    if not args.figures_only:
        for source in sources:
            start = time.monotonic()
            print(f'Executing {source.stem}', flush=True)
            subprocess.run([sys.executable, str(source)], cwd=ROOT, check=True)
            records.append({'section': source.stem, 'seconds': round(time.monotonic()-start, 2)})
    if args.figures or args.figures_only:
        subprocess.run([args.rscript, str(ROOT/'scripts/figures/render_all.R')], cwd=ROOT, check=True)
    stats = ROOT/'output/stats'
    stats.mkdir(parents=True, exist_ok=True)
    (stats/'last_run.json').write_text(json.dumps({
        'python': sys.version.split()[0],
        'packages': {p: importlib.metadata.version(p) for p in ('numpy', 'pandas', 'scipy', 'PyYAML')},
        'sections': records, 'figures': bool(args.figures or args.figures_only),
    }, indent=2)+'\n')


if __name__ == '__main__':
    main()
