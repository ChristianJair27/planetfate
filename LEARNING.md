# What we learned building planetfate (2026-11-09)

Notes by a self-taught team (a hobbyist and an AI agent), written so the next person can skip our detours.

## The physics, in the order it clicked

1. **Mass loss widens orbits, and the rule is embarrassingly simple.** If the star loses mass slowly compared with
   the orbital period, the planet's specific angular momentum sqrt(G M a (1-e^2)) is conserved while M drops, so
   a M = const: a_f = a_i M_i/M_f and e does not change. "Slowly" is the adiabatic condition Psi = |dM/dt|/(M n) << 1
   (Veras et al. 2011); for giant-branch winds (1e5-1e6 yr) and orbits of years it holds everywhere except for the
   widest orbits (hundreds of AU), where the planet can be lost instead. Our REBOUND run of the eight planets lands
   within 2 % of the rule.
2. **The initial-final mass relation is the exchange rate.** Cummings et al. (2018) measured it with white dwarfs in
   open clusters of known age: the Sun keeps 0.57 Msun (expansion x1.76), a 2 Msun star keeps 0.65 (x3.1), a 5 Msun
   star keeps 1.0 (x5). The kink at 2.85-3.6 Msun is the onset of helium ignition in a non-degenerate core.
3. **The surface is not the danger zone; the tide is.** A planet just outside a red giant raises a bulge on a star
   that rotates far slower than the planet orbits; the bulge lags and drains the orbit. Mustill & Villaver (2012)
   follow this through every AGB thermal pulse and find Jupiters pulled in from 3 AU (1 Msun) to 5 AU (5 Msun);
   lighter planets raise weaker tides and survive from 1.5-3 AU. We encode their anchors as a log-log interpolation,
   nothing more; the real curves depend on the tidal prescription and disagree at the 30 % level between papers.
4. **Put 2 and 3 together and you get a gap.** The closest a surviving planet can end up is a_crit M_i/M_f: 5 AU for
   the Sun's descendants, 25 AU for a 5 Msun progenitor. Nordhaus & Spiegel (2013) drew the same gap edge at 4-40 AU.
   Everything observed at 0.01-0.1 AU around white dwarfs is therefore a migrant.
5. **Why multi-planet systems go unstable after the star dies.** Hill stability depends on the spacing in mutual
   Hill radii, R_H ∝ a (m_p/M*)^(1/3). When the star loses mass every a grows by M_i/M_f but R_H grows by
   (M_i/M_f)^(4/3), so the spacing in Hill units shrinks by (M_f/M_i)^(1/3): from 4.6 to 3.8 for a Sun-like star
   (Debes & Sigurdsson 2002). Our runs: at Delta = 4.6 the control system (constant mass) stays circular for 2.5 Myr
   (e < 0.04); with mass loss the first close encounter comes at 0.46 Myr, before the mass loss is even finished,
   and one planet is thrown to e = 0.97 with a pericentre of 0.14 AU from the white dwarf. At Delta = 5.5 both runs
   stay stable over 2 Myr: the effect needs a system that was already near the edge, which is the paper's point.
6. **Why this matters for white-dwarf pollution.** A planet or planetesimal on an orbit grazing the white dwarf is
   tidally shredded (Roche limit ~0.005 AU) and the rock rains onto the star; that is the accepted origin of the
   metals in a quarter of all white dwarfs. The instability after mass loss is the delivery mechanism.

## Methods we copied from professionals

- IAS15 (Rein & Spiegel 2015) for everything with a varying central mass and close encounters: adaptive, symplectic
  enough for Myr runs, and it does not blow up when two Jupiters pass within a Hill radius. WHFast would be ten times
  faster for the quiet phases but is unsafe with M*(t) and encounters.
- Change the central mass between integrator steps (every 10 yr) rather than inside the force routine; with ramps of
  2e5 yr the result is indistinguishable from a continuous loss.
- Always run the control: the same system without mass loss, for the same total time. Without it the instability
  could be blamed on the initial conditions.
- Anchor the model on published numbers and name the paper for each; where papers disagree, say so and pick one.

## What we would do next

- Replace the four-anchor a_crit with a proper tidal integration (Zahn prescription, stellar tracks from MIST) so the
  limit depends continuously on planet mass, metallicity and the AGB history.
- Propagate the Archive's host-mass errors and the IFMR scatter by Monte Carlo per planet.
- Add the Kozai channel (a wide stellar companion) and tidal circularisation to follow a flung planet to 0.02 AU.
- Join with wdclock: for each of the 16,676 Gaia white dwarfs, the progenitor mass (inverse IFMR), hence the gap
  edge and the migrant-only region where TESS could see a transit.
