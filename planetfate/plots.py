"""Figures: dark style, provenance line, English or Spanish labels.

    F = Figures(lang="en", outdir="figures")
    F.ifmr()                      # initial-final mass relation and the expansion factor
    F.gap(known_wd_planets)       # a_crit and the closest surviving orbit vs progenitor mass; the "gap" to 0.01 AU
    F.known_planets(fates)        # the NASA Exoplanet Archive: who survives
    F.final_orbits(fates)         # where the survivors end up, against the white-dwarf habitable zone
    F.instability(run)            # REBOUND: the Debes & Sigurdsson effect
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

from .ifmr import final_mass
from .survival import a_crit

BG = "#070b14"
T = {
 "en": dict(prov="planetfate · IFMR: Cummings+ 2018 · engulfment limits: Mustill & Villaver 2012, Schröder & Connon Smith 2008 · planets: NASA Exoplanet Archive",
            ifmr_t="What is left of a star: the initial-final mass relation", mi="Initial stellar mass (Msun)", mf="White-dwarf mass (Msun)",
            exp="Orbit expansion factor M_i / M_f", gap_t="The gap: nothing survives the giant branches inside the red line, nothing ends up inside the orange one",
            a="Orbital distance (AU)", crit="engulfed or dragged in (Jupiter mass)", crit_e="engulfed (Earth mass)", amin="closest surviving orbit around the white dwarf",
            hz="white-dwarf habitable zone (0.005–0.02 AU)", known="known planets around white dwarfs", kp_t="The fate of {n:,} known exoplanets when their stars die",
            host="Host star mass (Msun)", surv="survives ({n})", eng="engulfed ({n})", alive="star outlives the Universe so far ({n})", sn="supernova ({n})",
            fo_t="Where the {n} survivors end up: orbits around the white dwarf", af="Final orbital distance (AU)", cnt="Number of planets",
            inst_t="Three Jupiters, stable for 300,000 years, then the star loses 43 % of its mass", tyr="Time (Myr)", ae="Semi-major axis (AU)", ecc="Eccentricity",
            loss="mass loss", enc="first close encounter"),
 "es": dict(prov="planetfate · IFMR: Cummings+ 2018 · límites de engullimiento: Mustill & Villaver 2012, Schröder & Connon Smith 2008 · planetas: NASA Exoplanet Archive",
            ifmr_t="Lo que queda de una estrella: la relación masa inicial–masa final", mi="Masa inicial de la estrella (Msol)", mf="Masa de la enana blanca (Msol)",
            exp="Factor de expansión de las órbitas M_i / M_f", gap_t="El hueco: nada sobrevive a las gigantes dentro de la línea roja; nada acaba dentro de la naranja",
            a="Distancia orbital (UA)", crit="engullido o arrastrado (masa de Júpiter)", crit_e="engullido (masa de la Tierra)", amin="órbita superviviente más cercana a la enana blanca",
            hz="zona habitable de la enana blanca (0,005–0,02 UA)", known="planetas conocidos alrededor de enanas blancas", kp_t="El destino de {n:,} exoplanetas conocidos cuando mueran sus estrellas",
            host="Masa de la estrella (Msol)", surv="sobrevive ({n})", eng="engullido ({n})", alive="la estrella durará más que el Universo hasta hoy ({n})", sn="supernova ({n})",
            fo_t="Dónde acaban los {n} supervivientes: órbitas alrededor de la enana blanca", af="Distancia orbital final (UA)", cnt="Número de planetas",
            inst_t="Tres Júpiter, estables 300.000 años; luego la estrella pierde el 43 % de su masa", tyr="Tiempo (Myr)", ae="Semieje mayor (UA)", ecc="Excentricidad",
            loss="pérdida de masa", enc="primer encuentro cercano"),
}
KNOWN_WD_PLANETS = [("WD 1856+534 b", 0.0204), ("WD J0914+1914 (giant debris)", 0.07), ("WD 0806-661 b", 2500.0), ("MOA-2010-BLG-477L b", 2.8), ("PSR B1620-26 b", 23.0)]


def style():
    plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.edgecolor": "#3a4256",
        "axes.labelcolor": "#dfe6f3", "xtick.color": "#aab4c8", "ytick.color": "#aab4c8", "text.color": "#eef2fa", "grid.color": "#1e2638",
        "axes.grid": True, "grid.linewidth": 0.5, "font.family": "DejaVu Sans", "axes.titleweight": "bold", "axes.titlesize": 13, "legend.frameon": False})


class Figures:
    def __init__(self, lang="en", outdir="figures"):
        self.t = T[lang]; self.lang = lang; self.out = Path(outdir); self.out.mkdir(parents=True, exist_ok=True); style()

    def _finish(self, fig, name, extra=""):
        fig.text(0.01, 0.005, self.t["prov"] + (" · " + extra if extra else ""), fontsize=7.5, color="#8a94a8")
        p = self.out / f"{name}.png"; fig.savefig(p, dpi=170, bbox_inches="tight"); plt.close(fig); return p

    def ifmr(self, name="ifmr"):
        m = np.linspace(0.8, 8, 400); mf = final_mass(m)
        fig, ax = plt.subplots(figsize=(9.5, 5.8)); ax2 = ax.twinx()
        ax.plot(m, mf, color="#9fd3ff", lw=2.5); ax.fill_between(m, mf - 0.06, mf + 0.06, color="#9fd3ff", alpha=0.15)
        ax2.plot(m, m / mf, color="#ffd27a", lw=2, ls="--"); ax2.set_ylabel(self.t["exp"], color="#ffd27a"); ax2.grid(False)
        for mi, lab in ((1.0, "Sun"), (2.0, "Sirius A"), (5.1, "Sirius B's progenitor")):
            ax.scatter([mi], final_mass(mi), color="white", zorder=5, s=30); ax.annotate(lab, (mi, final_mass(mi)[0]), xytext=(6, -14), textcoords="offset points", fontsize=9, color="#dfe6f3")
        ax.set_xlabel(self.t["mi"]); ax.set_ylabel(self.t["mf"]); ax.set_title(self.t["ifmr_t"]); ax.set_ylim(0.45, 1.35)
        return self._finish(fig, name, "band: 0.06 Msun scatter")

    def gap(self, name="gap"):
        m = np.linspace(0.8, 8, 300); mf = final_mass(m)
        cj = a_crit(m, 1.0); ce = a_crit(m, 1 / 317.8); amin = cj * m / mf
        fig, ax = plt.subplots(figsize=(10, 6.2))
        ax.fill_between(m, 1e-3, ce, color="#ff4d4d", alpha=0.25); ax.fill_between(m, ce, cj, color="#ff7a3d", alpha=0.18)
        ax.plot(m, cj, color="#ff6a3d", lw=2, label=self.t["crit"]); ax.plot(m, ce, color="#ff4d4d", lw=1.2, ls="--", label=self.t["crit_e"])
        ax.plot(m, amin, color="#ffd27a", lw=2.5, label=self.t["amin"])
        ax.axhspan(0.005, 0.02, color="#3fcf8e", alpha=0.35, label=self.t["hz"])
        for lab, a in KNOWN_WD_PLANETS[:2]:
            ax.scatter([0.9], [a], color="#4da3ff", s=40, zorder=6); ax.annotate(lab, (0.9, a), xytext=(8, 0), textcoords="offset points", fontsize=9, color="#9fd3ff", va="center")
        ax.set_yscale("log"); ax.set_ylim(3e-3, 60); ax.set_xlim(0.8, 8)
        ax.set_xlabel(self.t["mi"]); ax.set_ylabel(self.t["a"]); ax.set_title(self.t["gap_t"], fontsize=12); ax.legend(loc="upper left", fontsize=9)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        return self._finish(fig, name, "initial orbit below the red lines: engulfed; survivors end on the orange line or beyond (a_f = a_i M_i/M_f)")

    def known_planets(self, fates: pd.DataFrame, name="known_planets"):
        d = fates[fates.fate.isin(["survives", "engulfed", "star_still_alive", "supernova"]) & fates.pl_orbsmax.notna() & fates.st_mass.notna()]
        fig, ax = plt.subplots(figsize=(10.5, 6.5))
        cols = {"engulfed": ("#ff6a3d", self.t["eng"]), "star_still_alive": ("#5b6b86", self.t["alive"]), "survives": ("#3fcf8e", self.t["surv"]), "supernova": ("#ffffff", self.t["sn"])}
        for k in ("star_still_alive", "engulfed", "supernova", "survives"):
            s = d[d.fate == k]
            ax.scatter(s.st_mass, s.pl_orbsmax, s=9 if k != "survives" else 22, color=cols[k][0], alpha=0.75 if k != "engulfed" else 0.45, linewidths=0, label=cols[k][1].format(n=len(s)))
        m = np.linspace(0.8, 8, 200); ax.plot(m, a_crit(m, 1.0), color="#ff6a3d", lw=1.5, ls="--")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.07, 9); ax.set_ylim(0.005, 3000)
        ax.set_xlabel(self.t["host"]); ax.set_ylabel(self.t["a"]); ax.set_title(self.t["kp_t"].format(n=len(d))); ax.legend(loc="lower right", fontsize=9, markerscale=1.5)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        return self._finish(fig, name, "dashed: Jupiter-mass engulfment limit")

    def final_orbits(self, fates: pd.DataFrame, name="final_orbits"):
        s = fates[fates.fate == "survives"]
        fig, ax = plt.subplots(figsize=(10, 5.8))
        bins = np.logspace(-2.5, 4, 40)
        ax.hist(s.a_final_au, bins=bins, color="#3fcf8e", alpha=0.85)
        ax.axvspan(0.005, 0.02, color="#9fd3ff", alpha=0.35); ax.text(0.006, ax.get_ylim()[1] * 0.85, self.t["hz"], color="#9fd3ff", fontsize=9, rotation=90, va="top")
        for lab, a in KNOWN_WD_PLANETS[:1]:
            ax.axvline(a, color="#4da3ff", lw=1.5); ax.text(a * 1.15, ax.get_ylim()[1] * 0.6, lab, color="#9fd3ff", fontsize=9, rotation=90, va="top")
        ax.set_xscale("log"); ax.set_xlabel(self.t["af"]); ax.set_ylabel(self.t["cnt"]); ax.set_title(self.t["fo_t"].format(n=len(s)))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        return self._finish(fig, name, f"closest survivor: {s.a_final_au.min():.1f} AU; median {s.a_final_au.median():.0f} AU")

    def instability(self, run: dict, name="instability"):
        t = run["t"] / 1e6; a = run["a"]; e = run["e"]
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        cols = ["#9fd3ff", "#ffd27a", "#ff7a3d"]
        for i in range(a.shape[1]):
            ok = (e[:, i] < 1) & (a[:, i] > 0)
            ax1.plot(t[ok], a[ok, i], color=cols[i], lw=1.2); ax2.plot(t[ok], e[ok, i], color=cols[i], lw=1.2)
        t0, t1 = run["t_ms"] / 1e6, (run["t_ms"] + run["t_loss"]) / 1e6
        for ax in (ax1, ax2):
            ax.axvspan(t0, t1, color="#ff4d4d", alpha=0.15)
        ax1.text(t0, ax1.get_ylim()[1] * 0.92, self.t["loss"] + f" {run['m_star']:.2f} → {run['m_final']:.3f} Msun", color="#ff9a6a", fontsize=10)
        if run.get("first_encounter_yr"):
            for ax in (ax1, ax2):
                ax.axvline(run["first_encounter_yr"] / 1e6, color="white", lw=1, ls=":")
            ax2.text(run["first_encounter_yr"] / 1e6, 0.9, " " + self.t["enc"], color="white", fontsize=9)
        ax1.set_yscale("log"); ax1.set_ylabel(self.t["ae"]); ax2.set_ylabel(self.t["ecc"]); ax2.set_xlabel(self.t["tyr"]); ax2.set_ylim(0, 1)
        ax1.set_title(self.t["inst_t"])
        return self._finish(fig, name, f"REBOUND IAS15 · Δ = {run['delta']} mutual Hill radii · Debes & Sigurdsson 2002 effect")
