"""Command line.

    planetfate one --mstar 1.0 --a 5.2 --mp 1.0 [--e 0.05]      # fate of one planet (mp in Jupiter masses)
    planetfate ifmr --mstar 1 2 3 5                              # white-dwarf masses and expansion factors
    planetfate catalog exoplanets.csv fates.csv                  # NASA Exoplanet Archive pscomppars columns
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .ifmr import final_mass
from .survival import a_crit, fate, fate_table


def main(argv=None):
    ap = argparse.ArgumentParser(prog="planetfate", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    o = sub.add_parser("one"); o.add_argument("--mstar", type=float, required=True); o.add_argument("--a", type=float, required=True)
    o.add_argument("--mp", type=float, default=1.0, help="planet mass in Jupiter masses (Earth = 0.00315)"); o.add_argument("--e", type=float, default=0.0)
    i = sub.add_parser("ifmr"); i.add_argument("--mstar", type=float, nargs="+", default=[1, 1.5, 2, 3, 5])
    c = sub.add_parser("catalog"); c.add_argument("inp"); c.add_argument("out")
    a = ap.parse_args(argv)
    if a.cmd == "one":
        f = fate(a.mstar, a.a, a.mp, a.e)
        print(json.dumps({k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in f.__dict__.items()}, indent=1))
    elif a.cmd == "ifmr":
        for m in a.mstar:
            mf = final_mass(m)[0]
            print(f"M_i = {m:4.2f} Msun -> M_wd = {mf:.3f} Msun, orbits expand x{m / mf:.2f}, a_crit(Jupiter) = {a_crit(m, 1.0)[0]:.2f} AU, "
                  f"closest surviving Jupiter ends at {a_crit(m, 1.0)[0] * m / mf:.1f} AU")
    elif a.cmd == "catalog":
        import pandas as pd
        r = fate_table(pd.read_csv(a.inp)); r.to_csv(a.out, index=False)
        print(r["fate"].value_counts().to_string())


if __name__ == "__main__":
    sys.exit(main())
