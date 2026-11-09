"""planetfate: which planets survive the death of their star, and where they end up.

    from planetfate import fate, final_mass, a_crit
    fate(1.0, 5.2, m_planet_mjup=1.0)        # Jupiter: survives, a_f = 5.2 * 1.0/0.569 = 9.1 AU
    fate(1.0, 1.0, m_planet_mjup=1/317.8)    # Earth: engulfed (inside a_crit ~ 1.5 AU)
"""
from .ifmr import final_mass, initial_mass, mass_loss_fraction, SCATTER
from .survival import Fate, a_crit, fate, fate_table, period_years, r_max

__version__ = "0.1.0"
__all__ = ["final_mass", "initial_mass", "mass_loss_fraction", "SCATTER", "Fate", "a_crit", "fate", "fate_table", "period_years", "r_max"]
