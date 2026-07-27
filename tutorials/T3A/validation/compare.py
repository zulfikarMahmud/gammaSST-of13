#!/usr/bin/env python3
"""
Quantitative T3A validation for gammaSST.

gnuplot (./createGraphs) gives you the picture; this gives you the numbers,
so a regression can be caught without eyeballing a plot.

Reports:
  * transition onset location (Re_x where c_f first rises off the laminar curve)
  * c_f error vs the Savill/ERCOFTAC experimental data
  * freestream turbulence-intensity decay error
  * sanity bounds on gammaInt (must stay in [0, 1])

Usage:  python3 compare.py [caseDir]      (default: parent directory)
"""
import os
import re
import sys

import numpy as np

case = sys.argv[1] if len(sys.argv) > 1 else ".."

# --- T3A conditions -----------------------------------------------------------
UINF = 5.4        # m/s
NU = 1.5e-5       # m^2/s
XLE = 0.04        # leading edge location [m]


def latest_time(sub):
    d = os.path.join(case, "postProcessing", sub)
    if not os.path.isdir(d):
        sys.exit(f"[compare] missing {d} -- run the case first")
    times = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t)]
    return os.path.join(d, max(times, key=float), "line.xy")


# --- skin friction ------------------------------------------------------------
d = np.loadtxt(latest_time("wallShearStressGraph"))
x, tau_x = d[:, 0], d[:, 1]
Rex = (x - XLE) * UINF / NU
cf = -tau_x / (0.5 * UINF**2)

# --- experiment ---------------------------------------------------------------
e = np.loadtxt(os.path.join(os.path.dirname(__file__), "exptData", "T3A.dat"))
ex, ecf, etu = e[:, 0] / 1000.0, e[:, 1], e[:, 2]
eRex = ex * UINF / 1.51e-5

# --- transition onset ---------------------------------------------------------
# Standard proxy: the LOCAL c_f minimum, i.e. the end of the laminar decay and
# the start of the transitional rise.  The search window excludes the
# leading-edge region (where the c_f is singular / mesh dependent) and the far
# downstream turbulent region (where c_f decays again and would otherwise give
# a spurious global minimum).
valid = Rex > 1e4
onset_win = (Rex > 1e4) & (Rex < 3e5)
Rex_onset = Rex[onset_win][np.argmin(cf[onset_win])]

# experiment: same c_f-minimum proxy
eRex_onset = eRex[np.argmin(ecf)]

# --- errors on the common Re_x range -----------------------------------------
lo, hi = max(Rex.min(), eRex.min()), min(Rex.max(), eRex.max())
m = (Rex >= lo) & (Rex <= hi)
cf_i = np.interp(Rex[m], eRex, ecf)
err = cf[m] - cf_i
mae = np.mean(np.abs(err))
rel = np.mean(np.abs(err) / cf_i) * 100

# --- turbulence intensity decay ----------------------------------------------
dk = np.loadtxt(latest_time("kGraph"))
xk, kk = dk[:, 0], dk[:, 1]
tu = np.sqrt(2.0 / 3.0 * kk) / UINF * 100
tu_i = np.interp(xk - XLE, ex, etu)
tu_mae = np.mean(np.abs(tu - tu_i))

# --- gammaInt bounds ----------------------------------------------------------
def field_stats(name):
    times = [t for t in os.listdir(case) if re.fullmatch(r"[\d.]+", t)]
    if not times:
        return None
    p = os.path.join(case, max(times, key=float), name)
    if not os.path.isfile(p):
        return None
    txt = open(p).read()
    mm = re.search(r"internalField\s+nonuniform[^(]*\(\s*(.*?)\)\s*;", txt, re.S)
    if not mm:
        mm = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)", txt)
        v = float(mm.group(1))
        return v, v, v
    a = np.fromstring(mm.group(1), sep=" ")
    return a.min(), a.max(), a.mean()


print("=" * 66)
print("  T3A flat-plate transition -- gammaSST vs Savill/ERCOFTAC experiment")
print("=" * 66)
print(
    f"  transition onset  Re_x : {Rex_onset:.3e}   "
    f"(exp {eRex_onset:.3e}, ratio {Rex_onset/eRex_onset:.2f})"
)
print(f"  c_f  mean abs error    : {mae:.5f}   ({rel:.1f} % of experiment)")
print(f"  c_f  range   CFD       : {cf[valid].min():.5f} .. {cf[valid].max():.5f}")
print(f"  c_f  range   exp       : {ecf.min():.5f} .. {ecf.max():.5f}")
print(f"  Tu   mean abs error    : {tu_mae:.3f} %")

g = field_stats("gammaInt")
if g:
    lo_g, hi_g, mean_g = g
    # gamma is an intermittency: it must stay in [0, 1].  A sub-1% overshoot at
    # the discrete level is normal (the equation self-limits at 1 rather than
    # being clipped there); anything larger means the source-term signs or the
    # implicit/explicit split are wrong.
    ok = "OK" if (lo_g >= -1e-6 and hi_g <= 1.01) else "*** OUT OF BOUNDS ***"
    print(f"  gammaInt  min/max/mean : {lo_g:.4f} / {hi_g:.4f} / {mean_g:.4f}   {ok}")

k = field_stats("k")
if k:
    print(f"  k         min/max/mean : {k[0]:.4g} / {k[1]:.4g} / {k[2]:.4g}")
print("=" * 66)
