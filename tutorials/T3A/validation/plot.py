#!/usr/bin/env python3
"""
T3A validation figure: skin friction and freestream turbulence decay.

Plots the gammaSST result against the Savill/ERCOFTAC experimental data, and
optionally against a second case (e.g. the built-in kOmegaSSTLM) for comparison.

Usage:
    python3 plot.py                          # gammaSST vs experiment
    python3 plot.py --compare <caseDir>      # ... plus a second model
    python3 plot.py --compare <dir> --label kOmegaSSTLM
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

UINF = 5.4
NU = 1.5e-5
XLE = 0.04

# Categorical series colours (validated: CVD dE 24.7, normal dE 33.6).
C1 = "#2a78d6"   # slot 1 - primary model
C2 = "#eb6834"   # slot 2 - comparison model
INK = "#0b0b0b"          # experiment: ground truth drawn in primary ink
INK_2 = "#52514e"        # secondary ink for axes/labels
GRID = "#d9d8d4"


def latest(caseDir, sub):
    d = os.path.join(caseDir, "postProcessing", sub)
    if not os.path.isdir(d):
        sys.exit(f"[plot] missing {d} -- run the case first")
    times = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t)]
    return np.loadtxt(os.path.join(d, max(times, key=float), "line.xy"))


def cf_curve(caseDir):
    d = latest(caseDir, "wallShearStressGraph")
    Rex = (d[:, 0] - XLE) * UINF / NU
    cf = -d[:, 1] / (0.5 * UINF**2)
    m = Rex > 1e4
    return Rex[m], cf[m]


def tu_curve(caseDir):
    d = latest(caseDir, "kGraph")
    x = d[:, 0] - XLE
    tu = np.sqrt(2.0 / 3.0 * d[:, 1]) / UINF * 100
    return x, tu


ap = argparse.ArgumentParser()
ap.add_argument("case", nargs="?", default=os.path.join(HERE, ".."))
ap.add_argument("--compare", default=None)
ap.add_argument("--label", default="kOmegaSSTLM")
ap.add_argument("--out", default=os.path.join(HERE, "T3A_validation.png"))
a = ap.parse_args()

e = np.loadtxt(os.path.join(HERE, "exptData", "T3A.dat"))
eRex = e[:, 0] / 1000 * UINF / 1.51e-5
ecf = e[:, 1]
ex = e[:, 0] / 1000
etu = e[:, 2]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
fig.patch.set_facecolor("#fcfcfb")

for ax in (ax1, ax2):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9)

# ---- skin friction -----------------------------------------------------------
Rex, cf = cf_curve(a.case)
ax1.plot(Rex, cf, "-", color=C1, lw=2, label="gammaSST", zorder=3)
if a.compare:
    Rex2, cf2 = cf_curve(a.compare)
    ax1.plot(Rex2, cf2, "-", color=C2, lw=2, label=a.label, zorder=2)
ax1.plot(
    eRex, ecf, "o", mfc="none", mec=INK, mew=1.4, ms=6,
    label="Experiment (T3A)", zorder=4,
)
ax1.set_xlim(0, 5.5e5)
ax1.set_ylim(0, 0.008)
ax1.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
ax1.xaxis.get_offset_text().set_color(INK_2)
ax1.set_xlabel("$Re_x$", color=INK_2, fontsize=10)
ax1.set_ylabel("$c_f$", color=INK_2, fontsize=10)
ax1.set_title("Skin friction", color=INK, fontsize=11, loc="left", pad=8)
ax1.legend(frameon=False, fontsize=9, labelcolor=INK_2)

# ---- turbulence intensity decay ---------------------------------------------
x, tu = tu_curve(a.case)
ax2.plot(x, tu, "-", color=C1, lw=2, label="gammaSST", zorder=3)
if a.compare:
    x2, tu2 = tu_curve(a.compare)
    ax2.plot(x2, tu2, "-", color=C2, lw=2, label=a.label, zorder=2)
ax2.plot(
    ex, etu, "o", mfc="none", mec=INK, mew=1.4, ms=6,
    label="Experiment (T3A)", zorder=4,
)
ax2.set_xlim(0, 1.5)
ax2.set_ylim(0, 3.5)
ax2.set_xlabel("$x$  [m]", color=INK_2, fontsize=10)
ax2.set_ylabel("$Tu$  [%]", color=INK_2, fontsize=10)
ax2.set_title("Freestream turbulence decay", color=INK, fontsize=11, loc="left", pad=8)
ax2.legend(frameon=False, fontsize=9, labelcolor=INK_2)

fig.suptitle(
    "ERCOFTAC T3A flat plate  ·  $U_\\infty$ = 5.4 m/s, $Tu_{inlet}$ = 3.3 %",
    color=INK, fontsize=12, x=0.005, ha="left", y=0.99,
)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(a.out, dpi=140, facecolor=fig.get_facecolor())
print(f"[plot] wrote {a.out}")
