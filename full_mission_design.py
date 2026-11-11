"""Run the complete TARA_LV conceptual launch-to-docking design from subsystem CSVs."""
from pathlib import Path
import argparse
from tara.pipeline import run


def main():
    root=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,default=root/'inputs')
    p.add_argument('--output',type=Path,default=root/'results')
    p.add_argument('--profile',choices=['quick','standard'],default='standard',help='quick reduces Euler mesh/iterations and ascent targeting budget')
    p.add_argument('--no-animations',action='store_true');p.add_argument('--skip-cfd',action='store_true')
    p.add_argument('--skip-grid-study',action='store_true');p.add_argument('--reuse-ascent',action='store_true',help='reuse only when all CSV hashes match')
    p.add_argument('--strict',action='store_true',help='exit 2 if any design gate remains open')
    a=p.parse_args()
    summary=run(a.inputs,a.output,a.profile,not a.no_animations,not a.skip_cfd,not a.skip_grid_study,a.reuse_ascent)
    return 2 if a.strict and any(r['status']!='PASS' for r in summary['checks']) else 0


if __name__=='__main__':raise SystemExit(main())
