#!/usr/bin/env python3
"""
E387 Cp validation: gammaSST (OpenFOAM 13) vs McGhee NASA TM-4062 experiment.

pMean is computed HERE, in post, by averaging the per-write-time surface
samples written by the `surfaces` function object:

    postProcessing/surfaces/<time>/airfoilSurf.xy
        columns: x  y  z  p  tau_x  tau_y  tau_z

Why average at all?  The E387 laminar separation bubble is physically unsteady,
so a steady SIMPLE run limit-cycles instead of converging: on this case Cl
swings ~16 % peak-to-peak over the last 2000 iterations.  A single snapshot is
therefore not the answer -- the mean over the limit cycle is.

Cp = pMean / q  with  q = 0.5*Uinf^2  (OpenFOAM p is kinematic, p/rho).
Uinf = 3 -> q = 4.5.  Using the wrong q is the classic error here.

Because the average is built from a finite number of snapshots, the script also
reports the standard error of the mean and draws it as a band, so the sampling
uncertainty is visible rather than hidden.

Usage:  python3 validate_cp.py [caseDir] [Uinf] [tStart]
"""
import glob
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

case = sys.argv[1] if len(sys.argv) > 1 else "."
Uinf = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
tStart = float(sys.argv[3]) if len(sys.argv) > 3 else 2000.0
q = 0.5 * Uinf**2

expCsv = os.path.join(case, "validation", "exptData", "mcghee_e387_cp_alpha0.csv")

C_CFD = "#2a78d6"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

# ---- gather snapshots --------------------------------------------------------
files = sorted(
    glob.glob(os.path.join(case, "postProcessing/surfaces/*/airfoilSurf.xy")),
    key=lambda f: float(f.split(os.sep)[-2]),
)
files = [f for f in files if float(f.split(os.sep)[-2]) >= tStart]
if not files:
    sys.exit(f"[cp] no surface snapshots at t >= {tStart} under {case}")

times = [float(f.split(os.sep)[-2]) for f in files]
stack = [np.loadtxt(f) for f in files]
ref = stack[0]
for a in stack[1:]:                      # face ordering must be identical
    assert a.shape == ref.shape and np.allclose(a[:, 0], ref[:, 0])

P = np.array([a[:, 3] for a in stack])   # (nSnap, nFace)
pMean = P.mean(0)
pSEM = P.std(0, ddof=1) / np.sqrt(len(P))

x, y = ref[:, 0], ref[:, 1]
print(f"[cp] averaging {len(files)} snapshots  t = {times[0]:.0f} .. {times[-1]:.0f}   q = {q}")
print(f"[cp] snapshot-to-snapshot scatter: mean SEM = {pSEM.mean()/q:.4f} Cp,"
      f"  max = {pSEM.max()/q:.4f} Cp")


# ---- order points around the contour, split at the LE -----------------------
# E387 is cambered: a 'y > 0' split mixes surfaces near the LE, so walk the
# closed loop by nearest neighbour instead (same approach as lsb.py).
def walk(px, py):
    S = np.column_stack([px, py])
    order = [int(np.argmax(px))]
    un = set(range(len(S)))
    un.discard(order[0])
    while un:
        cur = S[order[-1]]
        idx = np.fromiter(un, int)
        j = idx[np.argmin(((S[idx] - cur) ** 2).sum(1))]
        order.append(int(j))
        un.discard(int(j))
    return np.array(order)


o = walk(x, y)
xs, ys, pm, pe = x[o], y[o], pMean[o], pSEM[o]
iLE = int(np.argmin(xs))
A, B = slice(0, iLE + 1), slice(iLE, len(xs))
upIsA = ys[A].mean() > ys[B].mean()

chord, x0 = x.max() - x.min(), x.min()


def surf(sl):
    xx, pp, ee = xs[sl], pm[sl], pe[sl]
    s = np.argsort(xx)
    return (xx[s] - x0) / chord, pp[s] / q, ee[s] / q


xc_u, cp_u, se_u = surf(A if upIsA else B)
xc_l, cp_l, se_l = surf(B if upIsA else A)

# ---- experiment --------------------------------------------------------------
e = np.genfromtxt(expCsv, delimiter=",", names=True)
ex, ecp_u, ecp_l = e["xc"], e["Cp_upper"], e["Cp_lower"]

mu = (ex >= xc_u.min()) & (ex <= xc_u.max())
ml = (ex >= xc_l.min()) & (ex <= xc_l.max())
du = np.interp(ex[mu], xc_u, cp_u) - ecp_u[mu]
dl = np.interp(ex[ml], xc_l, cp_l) - ecp_l[ml]

print(f"[cp] upper : mean|dCp| = {np.abs(du).mean():.4f}   max = {np.abs(du).max():.4f}")
print(f"[cp] lower : mean|dCp| = {np.abs(dl).mean():.4f}   max = {np.abs(dl).max():.4f}")
print(f"[cp] ALL   : mean|dCp| = {np.abs(np.r_[du, dl]).mean():.4f}")
print(f"[cp] suction peak: CFD {cp_u.min():.3f} @ x/c {xc_u[np.argmin(cp_u)]:.3f}"
      f"  |  exp {ecp_u.min():.3f} @ x/c {ex[np.argmin(ecp_u)]:.3f}")

# ---- plot --------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.8, 5.4))
fig.patch.set_facecolor("#fcfcfb")
ax.set_facecolor("#fcfcfb")
ax.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)

for xc, cp, se in ((xc_u, cp_u, se_u), (xc_l, cp_l, se_l)):
    ax.fill_between(xc, cp - 2 * se, cp + 2 * se, color=C_CFD, alpha=0.25, lw=0, zorder=2)
ax.plot(xc_u, cp_u, "-", color=C_CFD, lw=2, label="gammaSST (pMean, ±2 SEM)", zorder=3)
ax.plot(xc_l, cp_l, "-", color=C_CFD, lw=2, zorder=3)
ax.plot(ex, ecp_u, "o", mfc="none", mec=INK, mew=1.3, ms=5.5,
        label="McGhee NASA TM-4062", zorder=4)
ax.plot(ex, ecp_l, "o", mfc="none", mec=INK, mew=1.3, ms=5.5, zorder=4)

ax.invert_yaxis()
ax.set_xlabel("$x/c$", color=INK2, fontsize=10)
ax.set_ylabel("$C_p$", color=INK2, fontsize=10)
ax.set_title("E387  ·  Re = 200 000, $\\alpha$ = 0°, $Tu$ = 0.1 %  ·  mean $C_p$",
             color=INK, fontsize=11, loc="left", pad=8)
ax.legend(frameon=False, fontsize=9, labelcolor=INK2)
out = os.path.join(case, "validation", "cp_validation.png")
fig.tight_layout()
fig.savefig(out, dpi=140, facecolor=fig.get_facecolor())
print(f"[cp] plot -> {out}")
