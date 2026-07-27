#!/usr/bin/env python3
"""
Four-way T3A comparison at a single freestream velocity:

    this gammaSST  vs  furstj/gammaSST  vs  built-in kOmegaSSTLM  vs  experiment

All CFD cases must have been run at the SAME U (default 5.18 m/s, furstj's
value) so that the comparison against his digitised curve is like-for-like.

Usage:
    python3 plot_compare.py --ours <case> --lm <case> \
        [--furstj reference/furstj_T3A_cf.csv] [--uinf 5.18] [--out cmp.png]
"""
import argparse
import os
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
NU = 1.5e-5
XLE = 0.04

# Validated categorical hues (CVD dE 24.7 / normal dE 33.6 on adjacent pairs).
C_OURS = "#2a78d6"   # slot 1  blue
C_LM = "#eb6834"     # slot 2  orange
C_FUR = "#1baf7a"    # slot 3  aqua
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#d9d8d4"


def cf_curve(case, uinf):
    sub = os.path.join(case, "postProcessing", "wallShearStressGraph")
    t = max((d for d in os.listdir(sub) if re.fullmatch(r"[\d.eE+-]+", d)), key=float)
    d = np.loadtxt(os.path.join(sub, t, "line.xy"))
    Rex = (d[:, 0] - XLE) * uinf / NU
    cf = -d[:, 1] / (0.5 * uinf**2)
    m = Rex > 1e4
    return Rex[m], cf[m]


def onset(Rex, cf):
    w = (Rex > 1e4) & (Rex < 3e5)
    return Rex[w][np.argmin(cf[w])]


ap = argparse.ArgumentParser()
ap.add_argument("--ours", required=True)
ap.add_argument("--lm", required=True)
ap.add_argument("--furstj", default=os.path.join(HERE, "reference", "furstj_T3A_cf.csv"))
ap.add_argument("--uinf", type=float, default=5.18)
ap.add_argument("--out", default=os.path.join(HERE, "T3A_fourway.png"))
a = ap.parse_args()

oR, oC = cf_curve(a.ours, a.uinf)
lR, lC = cf_curve(a.lm, a.uinf)
f = np.loadtxt(a.furstj, delimiter=",")
fR, fC = f[:, 0], f[:, 1]

e = np.loadtxt(os.path.join(HERE, "exptData", "T3A.dat"))
eR = e[:, 0] / 1000 * a.uinf / 1.51e-5
eC = e[:, 1]

fig, (ax, axz) = plt.subplots(1, 2, figsize=(12.5, 4.8))
fig.patch.set_facecolor("#fcfcfb")

for x in (ax, axz):
    x.set_facecolor("#fcfcfb")
    x.grid(True, color=GRID, lw=0.6, zorder=0)
    x.set_axisbelow(True)
    for s in ("top", "right"):
        x.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        x.spines[s].set_color(GRID)
    x.tick_params(colors=INK2, labelsize=9)
    x.plot(fR, fC, "-", color=C_FUR, lw=2.4, label="gammaSST — furstj (digitised)", zorder=2)
    x.plot(lR, lC, "-", color=C_LM, lw=2, label="kOmegaSSTLM — built-in", zorder=3)
    x.plot(oR, oC, "-", color=C_OURS, lw=2, label="gammaSST — this repo", zorder=4)
    x.plot(eR, eC, "o", mfc="none", mec=INK, mew=1.4, ms=6.5, label="Experiment (T3A)", zorder=5)
    x.set_xlabel("$Re_x$", color=INK2, fontsize=10)
    x.set_ylabel("$c_f$", color=INK2, fontsize=10)
    x.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    x.xaxis.get_offset_text().set_color(INK2)

ax.set_xlim(0, 6e5)
ax.set_ylim(0, 0.008)
ax.set_title("Full range", color=INK, fontsize=11, loc="left", pad=8)
ax.legend(frameon=False, fontsize=9, labelcolor=INK2)

axz.set_xlim(5e4, 3.2e5)
axz.set_ylim(0.0015, 0.005)
axz.set_title("Transition region (detail)", color=INK, fontsize=11, loc="left", pad=8)

# mark onsets in the detail panel
for R_, C_, col in ((oR, oC, C_OURS), (lR, lC, C_LM), (fR, fC, C_FUR)):
    axz.axvline(onset(R_, C_), color=col, ls=":", lw=1.4, zorder=1)
axz.axvline(eR[np.argmin(eC)], color=INK, ls=":", lw=1.4, zorder=1)

fig.suptitle(
    f"ERCOFTAC T3A  ·  all CFD at $U_\\infty$ = {a.uinf} m/s, $Tu$ = 3.3 %, "
    r"$\nu_t/\nu$ = 12  ·  identical mesh & numerics",
    color=INK, fontsize=12, x=0.005, ha="left", y=0.99,
)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(a.out, dpi=140, facecolor=fig.get_facecolor())
print(f"[plot] wrote {a.out}")

print(f"\n  transition onset Re_x  (all at U = {a.uinf} m/s)")
print(f"    this repo gammaSST : {onset(oR, oC):.3e}")
print(f"    furstj    gammaSST : {onset(fR, fC):.3e}")
print(f"    built-in  LM       : {onset(lR, lC):.3e}")
print(f"    experiment         : {eR[np.argmin(eC)]:.3e}")
