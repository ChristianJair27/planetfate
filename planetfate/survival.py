"""Does a planet survive the death of its star, and where does it end up?

Three effects, in the order they act (Veras 2016, Royal Society Open Science 3,
150571, for the review):

1. **Engulfment and tides on the giant branches.** The star swells to about 1 AU
   (1 Msun; RGB tip 256 Rsun = 1.2 AU, Schroeder & Connon Smith 2008) and to over
   5 AU for a 5 Msun star on the AGB (Mustill & Villaver 2012). A planet inside the
   stellar surface is swallowed; a planet somewhat outside is dragged in by the
   tide it raises on the slowly rotating, hugely extended envelope, the more so the
   heavier the planet. Mustill & Villaver (2012, ApJ 761, 121), with full AGB
   thermal pulses and Zahn tides, find Jovian planets engulfed from initial orbits
   out to ~3 AU (1 Msun) and ~5 AU (5 Msun), terrestrial planets out to ~1.5 AU
   (1 Msun) and ~3 AU (5 Msun); Earth itself is lost near the RGB tip (Schroeder &
   Connon Smith 2008; survival would need a present orbit of about 1.15 AU).
   Debes & Sigurdsson (2002) put the solar-type limit at ~5 AU. We encode these
   published anchors as a critical initial semi-major axis a_crit(M_i, m_p),
   interpolated log-linearly in stellar mass and in planet mass between the
   "terrestrial" (1 Earth) and "Jovian" (1 Jupiter) anchors, and treat the pericentre
   a(1-e) as the distance that matters.

2. **Adiabatic orbit expansion.** Mass loss on the giant branches is slow compared
   with the orbit (the adiabatic regime, Psi = |dM/dt| / (M n) << 1, Veras et al.
   2011), so angular momentum per unit reduced mass is conserved and
   a_f = a_i * M_i / M_f, with the eccentricity unchanged (Debes & Sigurdsson 2002,
   eq. in Sec. 2.1; Veras 2016). With the IFMR this is a factor 1.7-5.

3. **Hill stability after mass loss.** For a multi-planet system the Hill spacing
   criterion Delta_c ~ sqrt(8/3 (e1^2+e2^2) + 9 mu^(2/3)) depends on mu = m_p/M*,
   which grows as the star loses mass. Systems that were marginally stable on the
   main sequence become unstable around the white dwarf, on timescales that shrink by
   one to two orders of magnitude (Debes & Sigurdsson 2002). Planets then scatter,
   some inward: the leading explanation for the rocky debris that pollutes a quarter
   of all white dwarfs. ``nbody.py`` demonstrates this with REBOUND.

Consequence for planet hunters: the closest a surviving planet can end up is
a_min = a_crit * M_i / M_f, which is 4-40 AU (Nordhaus & Spiegel 2013): orbital
periods of decades around the white dwarf. Anything found at 0.01-0.1 AU (as WD
1856+534 b at 0.02 AU) had to get there afterwards, by scattering or Kozai cycles
plus tidal circularisation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ifmr import final_mass, M_MIN_WD_PROGENITOR, M_MAX_WD_PROGENITOR

MEARTH_IN_MJUP = 1.0 / 317.8
G_CGS = 6.674e-8
MSUN_G = 1.98847e33
AU_CM = 1.495978707e13
YEAR_S = 3.15576e7

# Critical initial semi-major axis for engulfment [AU], anchors from the literature (see module docstring):
# (stellar mass Msun, a_crit for a 1 Earth-mass planet, a_crit for a 1 Jupiter-mass planet)
ANCHORS = [
    (1.0, 1.5, 3.0),     # Mustill & Villaver 2012 (abstract): terrestrial ~1.5 AU, Jovian ~3 AU
    (5.0, 3.0, 5.0),     # Mustill & Villaver 2012 (abstract): terrestrial ~3 AU, Jovian ~5 AU
]
# Maximum stellar radii [AU] used for the "inside the star" flag (approximate): Schroeder & Connon Smith 2008 give
# 1.2 AU at the RGB tip for the Sun; Mustill & Villaver 2012 say a 5 Msun star "expands to over 5 AU" on the AGB.
R_MAX_ANCHORS = [(1.0, 1.2), (5.0, 5.0)]


def _loginterp(m, anchors_m, anchors_v):
    lm = np.log10(np.clip(m, 0.1, 20.0))
    return 10 ** np.interp(lm, np.log10(anchors_m), np.log10(anchors_v))


def a_crit(m_star, m_planet_mjup):
    """Critical initial semi-major axis [AU] below which the planet is engulfed or tidally dragged in.

    Interpolates between the Earth-mass and Jupiter-mass anchors in log(planet mass), and between the
    1 and 5 Msun anchors in log(stellar mass); outside those ranges the nearest anchor is held. Planets
    heavier than Jupiter use the Jupiter value (the dependence flattens; brown dwarfs are not handled).
    """
    m_star = np.atleast_1d(np.asarray(m_star, float)); mp = np.atleast_1d(np.asarray(m_planet_mjup, float))
    m_star, mp = np.broadcast_arrays(m_star, mp)
    am = np.array([a[0] for a in ANCHORS]); a_e = np.array([a[1] for a in ANCHORS]); a_j = np.array([a[2] for a in ANCHORS])
    ce = _loginterp(m_star, am, a_e); cj = _loginterp(m_star, am, a_j)
    w = np.clip((np.log10(np.clip(mp, 1e-4, None)) - np.log10(MEARTH_IN_MJUP)) / (0 - np.log10(MEARTH_IN_MJUP)), 0, 1)
    return ce * (1 - w) + cj * w


def r_max(m_star):
    """Approximate maximum giant radius [AU] (RGB/AGB) of the progenitor."""
    m = np.atleast_1d(np.asarray(m_star, float))
    return _loginterp(m, np.array([r[0] for r in R_MAX_ANCHORS]), np.array([r[1] for r in R_MAX_ANCHORS]))


def period_years(a_au, m_msun):
    a = np.asarray(a_au, float) * AU_CM
    return 2 * np.pi * np.sqrt(a ** 3 / (G_CGS * np.asarray(m_msun, float) * MSUN_G)) / YEAR_S


@dataclass
class Fate:
    verdict: str          # 'star_still_alive' | 'supernova' | 'engulfed' | 'survives'
    m_initial: float
    m_final: float
    a_initial: float
    a_final: float        # AU (NaN unless survives)
    period_final_yr: float
    a_crit: float
    expansion: float      # a_f / a_i


def fate(m_star, a_au, m_planet_mjup=1.0, e=0.0):
    """Fate of one planet: engulfed on the giant branches, or surviving on an expanded orbit."""
    m_star, a_au, mp, e = (float(x) for x in (m_star, a_au, m_planet_mjup, e))
    if m_star < M_MIN_WD_PROGENITOR:
        return Fate("star_still_alive", m_star, np.nan, a_au, np.nan, np.nan, np.nan, np.nan)
    if m_star > M_MAX_WD_PROGENITOR:
        return Fate("supernova", m_star, np.nan, a_au, np.nan, np.nan, np.nan, np.nan)
    mf = float(final_mass(m_star)[0]); ac = float(a_crit(m_star, mp)[0])
    peri = a_au * (1 - e)
    if peri < ac:
        return Fate("engulfed", m_star, mf, a_au, np.nan, np.nan, ac, m_star / mf)
    af = a_au * m_star / mf
    return Fate("survives", m_star, mf, a_au, af, float(period_years(af, mf)), ac, m_star / mf)


def fate_table(df, m_star="st_mass", a="pl_orbsmax", mp="pl_bmassj", e="pl_orbeccen"):
    """Vectorised fates for a DataFrame (NASA Exoplanet Archive column names by default)."""
    import pandas as pd
    out = df.copy()
    ms = out[m_star].to_numpy(float); aa = out[a].to_numpy(float)
    mpj = out[mp].to_numpy(float) if mp in out else np.full(len(out), 1.0)
    ee = np.nan_to_num(out[e].to_numpy(float), nan=0.0) if e in out else np.zeros(len(out))
    mf = final_mass(ms); ac = a_crit(ms, np.nan_to_num(mpj, nan=1.0))
    verdict = np.full(len(out), "no_data", dtype=object)
    ok = np.isfinite(ms) & np.isfinite(aa)
    verdict[ok & (ms < M_MIN_WD_PROGENITOR)] = "star_still_alive"
    verdict[ok & (ms > M_MAX_WD_PROGENITOR)] = "supernova"
    wd = ok & (ms >= M_MIN_WD_PROGENITOR) & (ms <= M_MAX_WD_PROGENITOR)
    peri = aa * (1 - ee)
    verdict[wd & (peri < ac)] = "engulfed"
    surv = wd & (peri >= ac)
    verdict[surv] = "survives"
    out["wd_mass"] = np.where(wd, mf, np.nan)
    out["a_crit_au"] = np.where(wd, ac, np.nan)
    out["a_final_au"] = np.where(surv, aa * ms / mf, np.nan)
    out["period_final_yr"] = np.where(surv, period_years(np.nan_to_num(out["a_final_au"], nan=1.0), np.nan_to_num(mf, nan=1.0)), np.nan)
    out["fate"] = verdict
    return out
