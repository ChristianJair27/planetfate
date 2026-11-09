"""Initial-final mass relation (IFMR): how much of a star is left as a white dwarf.

Physics
-------
A star below about 8 solar masses ends as a white dwarf made of the carbon-oxygen
core it built during helium burning; everything outside that core is blown away on
the red-giant and asymptotic-giant branches. The core mass grows with the initial
mass, slowly, so a 1 Msun star leaves 0.55-0.6 Msun while a 5 Msun star leaves about
1.0 Msun: the heavier the star, the larger the fraction lost, and (see
``survival.py``) the more its planets' orbits expand.

The relation used here is the semi-empirical, MIST-based fit of Cummings, Kalirai,
Tremblay, Ramirez-Ruiz & Choi (2018, ApJ 866, 21; arXiv:1809.01673), from white
dwarfs in open clusters of known age. Three linear segments (their eqs. 4-6):

    0.83 < M_i < 2.85 Msun :  M_f = 0.080 M_i + 0.489
    2.85 < M_i < 3.60 Msun :  M_f = 0.187 M_i + 0.184
    3.60 < M_i < 7.20 Msun :  M_f = 0.107 M_i + 0.471

with a scatter of 0.06 Msun. The segments meet at 0.717 Msun (2.85 Msun) and
0.857 Msun (3.60 Msun). Outside 0.83-7.20 Msun the relation is extrapolated and
flagged: below ~0.8 Msun the star has not had time to leave the main sequence in
the age of the Universe; above ~8 Msun it explodes as a supernova instead.
"""
from __future__ import annotations

import numpy as np

SEGMENTS = [(0.83, 2.85, 0.080, 0.489), (2.85, 3.60, 0.187, 0.184), (3.60, 7.20, 0.107, 0.471)]
SCATTER = 0.06          # Msun, Cummings et al. 2018
M_MIN_WD_PROGENITOR = 0.80   # below this the star is still on the main sequence after 13.8 Gyr (approx.)
M_MAX_WD_PROGENITOR = 8.0    # above this: core-collapse supernova (approx.)


def final_mass(m_initial):
    """White-dwarf mass [Msun] for initial mass(es) [Msun]; NaN outside 0.8-8 Msun."""
    m = np.atleast_1d(np.asarray(m_initial, float))
    out = np.full(m.shape, np.nan)
    for lo, hi, slope, icpt in SEGMENTS:
        sel = (m >= lo) & (m < hi)
        out[sel] = slope * m[sel] + icpt
    # gentle extrapolation of the end segments inside the physical window
    lo_sel = (m >= M_MIN_WD_PROGENITOR) & (m < SEGMENTS[0][0])
    out[lo_sel] = SEGMENTS[0][2] * m[lo_sel] + SEGMENTS[0][3]
    hi_sel = (m >= SEGMENTS[-1][1]) & (m <= M_MAX_WD_PROGENITOR)
    out[hi_sel] = SEGMENTS[-1][2] * m[hi_sel] + SEGMENTS[-1][3]
    return out


def mass_loss_fraction(m_initial):
    """Fraction of the star's mass lost before it becomes a white dwarf."""
    m = np.atleast_1d(np.asarray(m_initial, float))
    return 1.0 - final_mass(m) / m


def initial_mass(m_final):
    """Invert the relation (monotonic): progenitor mass for a white-dwarf mass; NaN outside 0.55-1.24 Msun."""
    mf = np.atleast_1d(np.asarray(m_final, float))
    grid = np.linspace(M_MIN_WD_PROGENITOR, M_MAX_WD_PROGENITOR, 2000)
    fg = final_mass(grid)
    out = np.interp(mf, fg, grid, left=np.nan, right=np.nan)
    return out
