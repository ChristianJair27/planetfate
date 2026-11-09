"""Checks against published numbers.

Cummings et al. 2018: the three IFMR segments meet at 0.717 Msun (2.85 Msun) and 0.857 Msun (3.60 Msun); a 1 Msun
star leaves ~0.57 Msun. Schroeder & Connon Smith 2008: the Sun ends at ~0.54 Msun (their model; the IFMR gives 0.57,
within the 0.06 scatter) and Earth at 1 AU is engulfed. Mustill & Villaver 2012: Jupiter-mass planets at 3 AU
(1 Msun) and 5 AU (5 Msun) are at the survival limit. Debes & Sigurdsson 2002 / Veras 2016: a_f = a_i M_i/M_f.
Nordhaus & Spiegel 2013: the closest surviving planets end 4-40 AU from the white dwarf.
"""
import numpy as np
import pytest

from planetfate import a_crit, fate, final_mass, initial_mass
from planetfate.nbody import hill_spacing_system, mass_loss_run


def test_ifmr_segments_meet():
    assert abs(final_mass(2.85)[0] - 0.717) < 0.005
    assert abs(final_mass(3.60)[0] - 0.857) < 0.005


def test_sun_final_mass_within_scatter():
    assert abs(final_mass(1.0)[0] - 0.546) < 0.06


def test_ifmr_inverse():
    assert abs(initial_mass(final_mass(2.0))[0] - 2.0) < 0.01


def test_earth_engulfed_jupiter_survives():
    assert fate(1.0, 1.0, 1 / 317.8).verdict == "engulfed"
    j = fate(1.0, 5.2, 1.0)
    assert j.verdict == "survives"
    assert abs(j.a_final - 5.2 * 1.0 / final_mass(1.0)[0]) < 1e-6


def test_acrit_anchors():
    assert abs(a_crit(1.0, 1.0)[0] - 3.0) < 1e-6 and abs(a_crit(5.0, 1.0)[0] - 5.0) < 1e-6
    assert abs(a_crit(1.0, 1 / 317.8)[0] - 1.5) < 1e-6


def test_minimum_wd_orbit_in_nordhaus_spiegel_range():
    for m in (1.0, 2.0, 3.0, 5.0):
        amin = a_crit(m, 1.0)[0] * m / final_mass(m)[0]
        assert 4.0 <= amin <= 40.0


def test_hill_spacing_monotonic():
    a = hill_spacing_system(1.0, 1e-3, 5.0, 5.5, 3)
    assert np.all(np.diff(a) > 0) and a[0] == 5.0


@pytest.mark.slow
def test_adiabatic_expansion_nbody():
    r = mass_loss_run([(1.0, 0.0, 3e-6)], m_initial=1.5, m_final=0.52, t_loss_yr=2e4, coast_yr=2e3, dt_update_yr=2.0, n_out=50)
    assert abs(r["a_final_nbody"][0] / r["a_final_adiabatic"][0] - 1) < 0.02
