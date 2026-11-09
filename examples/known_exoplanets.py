"""Worked example: the fate of every known exoplanet.

    python examples/known_exoplanets.py [--download] [--lang en|es] [--out figures]

Downloads (or reuses) the NASA Exoplanet Archive composite table (pscomppars) through its TAP service,
applies ``planetfate.fate_table`` and writes data/exoplanet_fates.csv plus the figures. The TAP query:

    select pl_name,hostname,discoverymethod,disc_year,sy_dist,st_mass,st_rad,st_teff,st_lum,pl_orbsmax,
           pl_orbper,pl_orbeccen,pl_bmasse,pl_bmassj,pl_rade,sy_snum,sy_pnum from pscomppars
"""
import argparse
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from planetfate import fate_table  # noqa: E402
from planetfate.plots import Figures  # noqa: E402

TAP = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,hostname,discoverymethod,disc_year,sy_dist,"
       "st_mass,st_rad,st_teff,st_lum,pl_orbsmax,pl_orbper,pl_orbeccen,pl_bmasse,pl_bmassj,pl_rade,sy_snum,sy_pnum+from+pscomppars&format=csv")

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", default="data/exoplanets_pscomppars.csv")
ap.add_argument("--download", action="store_true")
ap.add_argument("--lang", default="en"); ap.add_argument("--out", default="figures")
a = ap.parse_args()
cat = Path(a.catalog)
if a.download or not cat.exists():
    cat.parent.mkdir(parents=True, exist_ok=True); urllib.request.urlretrieve(TAP, cat); print("downloaded", cat)
df = pd.read_csv(cat)
r = fate_table(df)
r.to_csv(cat.with_name("exoplanet_fates.csv"), index=False)
counts = r["fate"].value_counts().to_dict()
s = r[r.fate == "survives"]
summary = {"planets": int(len(df)), **{k: int(v) for k, v in counts.items()},
           "engulfed_fraction_of_wd_hosts": round(float((r.fate == "engulfed").sum() / (r.fate.isin(["engulfed", "survives"])).sum()), 3),
           "closest_survivor_au": round(float(s.a_final_au.min()), 2), "median_survivor_au": round(float(s.a_final_au.median()), 1),
           "median_survivor_period_yr": round(float(s.period_final_yr.median()), 1),
           "survive_pct_by_method": {k: round(float(v) * 100, 1) for k, v in r[r.fate.isin(["engulfed", "survives"])].groupby("discoverymethod")["fate"].apply(lambda x: (x == "survives").mean()).items()}}
print(json.dumps(summary, indent=1))
F = Figures(lang=a.lang, outdir=a.out)
for fn in (F.ifmr, F.gap): print(fn())
print(F.known_planets(r)); print(F.final_orbits(r))
