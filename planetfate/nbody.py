"""Direct N-body checks with REBOUND (Rein & Liu 2012).

Two experiments:

``mass_loss_run``   one or more planets around a star whose mass is ramped down linearly
                    (the adiabatic limit: ramps of 1e5-1e6 yr against orbital periods of
                    years). Checks a * M = const and returns the orbital history.

``instability_run`` the Debes & Sigurdsson (2002) effect. Three equal-mass planets are
                    placed at a fixed mutual Hill spacing Delta (marginally stable on the
                    main sequence), integrated for a while at constant stellar mass, then
                    the star loses mass to its white-dwarf value. Because the mutual Hill
                    radius scales with (m_p/M*)^(1/3) and the spacing in Hill radii shrinks
                    when M* drops, the same system becomes unstable: eccentricities grow,
                    orbits cross, a planet is scattered inward or ejected. The run records
                    the first close encounter after mass loss and the minimum pericentre
                    reached, which is the quantity that matters for white-dwarf pollution.

Units: AU, yr, Msun (G = 4 pi^2). Integrator: IAS15 (adaptive, handles close encounters);
WHFast is faster for the constant-mass phases but cannot take the varying central mass
safely, so IAS15 throughout.
"""
from __future__ import annotations

import time

import numpy as np

from .ifmr import final_mass


def _sim():
    import rebound
    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "ias15"
    return sim


def mass_loss_run(planets, m_initial=1.0, m_final=None, t_loss_yr=2e5, coast_yr=2e4, dt_update_yr=5.0, n_out=400, seed=7, log=None):
    """planets: list of (a_au, e, m_msun). Returns dict with t, m, a[n_t, n_p], e[n_t, n_p], expansion check."""
    import rebound  # noqa: F401
    m_final = float(final_mass(m_initial)[0]) if m_final is None else m_final
    sim = _sim(); sim.add(m=m_initial)
    rng = np.random.default_rng(seed)
    for a, e, m in planets:
        sim.add(m=m, a=a, e=e, inc=np.radians(rng.uniform(0, 2)), Omega=rng.uniform(0, 2 * np.pi), omega=rng.uniform(0, 2 * np.pi), f=rng.uniform(0, 2 * np.pi))
    sim.move_to_com()
    total = t_loss_yr + coast_yr
    n_updates = int(round(total / dt_update_yr)); every = max(1, n_updates // n_out)
    rows = []; t0 = time.time()
    for k in range(n_updates + 1):
        t = k * dt_update_yr
        frac = min(1.0, t / t_loss_yr)
        sim.particles[0].m = m_initial + (m_final - m_initial) * frac
        if k % every == 0 or k == n_updates:
            orbs = sim.orbits(primary=sim.particles[0])
            rows.append([t, sim.particles[0].m] + [o.a for o in orbs] + [o.e for o in orbs])
            if log and k % (every * 20) == 0:
                log(f"  t={t:.0f} yr M={sim.particles[0].m:.3f} a={[round(o.a, 3) for o in orbs]}")
        if k < n_updates:
            sim.integrate(t + dt_update_yr)
    arr = np.array(rows); n = len(planets)
    a0 = np.array([p[0] for p in planets])
    return {"t": arr[:, 0], "m": arr[:, 1], "a": arr[:, 2:2 + n], "e": arr[:, 2 + n:2 + 2 * n],
            "a_final_nbody": arr[-1, 2:2 + n], "a_final_adiabatic": a0 * m_initial / m_final,
            "m_initial": m_initial, "m_final": m_final, "wall_s": time.time() - t0}


def hill_spacing_system(m_star=1.0, m_planet=1e-3, a_inner=5.0, delta=5.5, n_planets=3):
    """Semi-major axes of n equal-mass planets separated by `delta` mutual Hill radii (Chambers et al. 1996 convention)."""
    a = [a_inner]
    for _ in range(n_planets - 1):
        a1 = a[-1]
        # R_H,mutual = ((m1+m2)/(3 M))^(1/3) (a1+a2)/2 ; solve a2 - a1 = delta * R_H
        k = delta * ((2 * m_planet) / (3 * m_star)) ** (1 / 3) / 2
        a.append(a1 * (1 + k) / (1 - k))
    return np.array(a)


def instability_run(m_star=1.0, m_planet=1e-3, a_inner=5.0, delta=5.5, t_ms_yr=2e5, t_loss_yr=2e5, t_wd_yr=1.5e6,
                    dt_update_yr=10.0, n_out=1500, seed=11, log=None):
    """Three Jupiter-mass planets at Hill spacing `delta`: constant-mass phase, mass loss, white-dwarf phase.

    Returns the orbital history plus diagnostics: time of first close encounter (< 3 mutual Hill radii) after mass loss,
    minimum pericentre, planets ejected (e >= 1 or a > 1e4 AU).
    """
    import rebound  # noqa: F401
    m_final = float(final_mass(m_star)[0])
    a_list = hill_spacing_system(m_star, m_planet, a_inner, delta)
    sim = _sim(); sim.add(m=m_star)
    rng = np.random.default_rng(seed)
    for a in a_list:
        sim.add(m=m_planet, a=a, e=0.01, inc=np.radians(rng.uniform(0, 1)), Omega=rng.uniform(0, 2 * np.pi), omega=rng.uniform(0, 2 * np.pi), f=rng.uniform(0, 2 * np.pi))
    sim.move_to_com(); sim.exit_max_distance = 2e4
    total = t_ms_yr + t_loss_yr + t_wd_yr
    n_updates = int(round(total / dt_update_yr)); every = max(1, n_updates // n_out)
    rows = []; t0 = time.time(); first_enc = None; min_peri = np.inf; ejected = set()
    for k in range(n_updates + 1):
        t = k * dt_update_yr
        if t <= t_ms_yr:
            m_now = m_star
        elif t <= t_ms_yr + t_loss_yr:
            m_now = m_star + (m_final - m_star) * (t - t_ms_yr) / t_loss_yr
        else:
            m_now = m_final
        sim.particles[0].m = m_now
        orbs = sim.orbits(primary=sim.particles[0])
        for i, o in enumerate(orbs):
            if o.e < 1 and o.a > 0:
                min_peri = min(min_peri, o.a * (1 - o.e))
            if o.e >= 1 or o.a < 0 or o.a > 1e4:
                ejected.add(i)
        if first_enc is None and t > t_ms_yr:
            ps = sim.particles
            for i in range(1, sim.N):
                for j in range(i + 1, sim.N):
                    d = np.sqrt((ps[i].x - ps[j].x) ** 2 + (ps[i].y - ps[j].y) ** 2 + (ps[i].z - ps[j].z) ** 2)
                    rh = ((ps[i].m + ps[j].m) / (3 * m_now)) ** (1 / 3) * 0.5 * (abs(orbs[i - 1].a) + abs(orbs[j - 1].a))
                    if d < 3 * rh:
                        first_enc = t
        if k % every == 0 or k == n_updates:
            rows.append([t, m_now] + [o.a for o in orbs] + [o.e for o in orbs] + [o.inc for o in orbs])
            if log and k % (every * 50) == 0:
                log(f"  t={t / 1e6:.2f} Myr M={m_now:.3f} a={[round(o.a, 2) for o in orbs]} e={[round(o.e, 2) for o in orbs]}")
        if k < n_updates:
            try:
                sim.integrate(t + dt_update_yr)
            except Exception as ex:   # rebound.Escape
                if log:
                    log(f"  escape at t={t:.0f} yr: {ex}")
                break
    arr = np.array(rows); n = len(a_list)
    return {"t": arr[:, 0], "m": arr[:, 1], "a": arr[:, 2:2 + n], "e": arr[:, 2 + n:2 + 2 * n], "inc": arr[:, 2 + 2 * n:2 + 3 * n],
            "a_initial": a_list, "delta": delta, "m_star": m_star, "m_final": m_final, "t_ms": t_ms_yr, "t_loss": t_loss_yr,
            "first_encounter_yr": first_enc, "min_pericentre_au": float(min_peri), "ejected": sorted(ejected), "wall_s": time.time() - t0}
