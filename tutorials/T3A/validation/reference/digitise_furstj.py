#!/usr/bin/env python3
"""
Digitise the published c_f curve from furstj/gammaSST's T3A validation figure.

WHY: furstj's repository ships only figures for his T3A result, not the
underlying numbers, and his model targets the ESI OpenFOAM line so it cannot be
compiled against OpenFOAM 13 to reproduce them directly.  To compare against his
result quantitatively rather than by eye, we recover the curve from the figure.

The extraction is SELF-CALIBRATING: his figure also plots the experimental
points, and we hold that experimental data exactly (it is byte-identical to
OpenFOAM 13's own copy).  Predicting where those points must land and checking
them against the detected markers verifies the axis mapping before any curve
value is trusted -- it currently agrees to a mean of 0.44 px (max 0.82 px).

Source figure (not redistributed here -- clone his repo to regenerate):
    https://github.com/furstj/gammaSST
    testCases/T3A/validation/figures/Rex_vs_cf.png

Usage:
    python3 digitise_furstj.py /path/to/gammaSST/testCases/T3A/validation/figures/Rex_vs_cf.png
"""
import sys

import numpy as np

try:
    from PIL import Image
except ImportError:
    sys.exit("[digitise] needs Pillow:  pip install pillow")

OUT = "furstj_T3A_cf.csv"

# Axis ranges are set explicitly by his gnuplot script:
#     plot [:6e+5][0:0.01] ...
XMIN, XMAX = 0.0, 6.0e5      # Re_x
YMIN, YMAX = 0.0, 0.01       # c_f
UINF_FURSTJ = 5.18           # his freestream velocity
CURVE_RGB = (148, 0, 211)    # gnuplot pngcairo default linetype 1 (dark violet)


def main(path, expt):
    a = np.array(Image.open(path).convert("RGB")).astype(int)

    # --- locate the plot frame (long runs of dark pixels) ---------------------
    dark = a.max(2) < 100
    rows = np.nonzero(dark.sum(1) > 400)[0]
    cols = np.nonzero(dark.sum(0) > 400)[0]
    T, B = int(rows.min()), int(rows.max())
    L, R = int(cols.min()), int(cols.max())

    def px2x(px):
        return XMIN + (px - L) / (R - L) * (XMAX - XMIN)

    def px2y(py):
        return YMAX - (py - T) / (B - T) * (YMAX - YMIN)

    def x2px(x):
        return L + (x - XMIN) / (XMAX - XMIN) * (R - L)

    def y2py(y):
        return T + (YMAX - y) / (YMAX - YMIN) * (B - T)

    # --- calibration check against the known experimental markers -------------
    inside = np.zeros_like(dark)
    inside[T + 3:B - 2, L + 3:R - 2] = True
    inside[T:T + 120, R - 330:R] = False           # legend box
    ys, xs = np.nonzero(dark & inside)
    pts = np.column_stack([xs, ys])

    used = np.zeros(len(pts), bool)
    cents = []
    for i in range(len(pts)):
        if used[i]:
            continue
        grp = (np.abs(pts - pts[i]).max(1) < 9) & (~used)
        used |= grp
        cents.append(pts[grp].mean(0))
    cents = np.array(cents)

    e = np.loadtxt(expt)
    eRex = e[:, 0] / 1000 * UINF_FURSTJ / 1.51e-5
    pred = np.column_stack([x2px(eRex), y2py(e[:, 1])])
    pred = pred[(pred[:, 0] >= L) & (pred[:, 0] <= R)]

    err = [np.hypot(*(cents[np.argmin(np.abs(cents - p).sum(1))] - p)) for p in pred]
    err = np.array(err)
    print(f"[digitise] frame x[{L},{R}] y[{T},{B}]")
    print(f"[digitise] calibration on {len(pred)} known experimental markers:")
    print(f"           mean {err.mean():.2f} px, max {err.max():.2f} px")
    print(f"           1 px = {(XMAX-XMIN)/(R-L):.0f} Re_x, {(YMAX-YMIN)/(B-T):.2e} c_f")
    if err.mean() > 2.0:
        sys.exit("[digitise] calibration FAILED -- axis mapping is wrong")

    # --- extract the curve ----------------------------------------------------
    curve = np.abs(a - np.array(CURVE_RGB)).sum(2) < 180
    curve[:T + 1] = False
    curve[B:] = False
    curve[T:T + 40, R - 330:R] = False              # legend key sample

    rec = []
    for px in range(L + 1, R):
        yy = np.nonzero(curve[:, px])[0]
        if yy.size:
            rec.append((px2x(px), px2y(yy.mean())))
    rec = np.array(rec)

    hdr = (
        "Digitised from furstj/gammaSST T3A validation figure (Rex_vs_cf.png).\n"
        f"Freestream U = {UINF_FURSTJ} m/s, Tu = 3.3 %, nut/nu = 12.\n"
        f"Calibrated on the experimental markers: mean {err.mean():.2f} px "
        f"(1 px ~ {(XMAX-XMIN)/(R-L):.0f} in Re_x, {(YMAX-YMIN)/(B-T):.2e} in c_f).\n"
        "Re_x,cf"
    )
    np.savetxt(OUT, rec, delimiter=",", header=hdr, fmt="%.6g")
    print(f"[digitise] wrote {OUT}  ({len(rec)} points, "
          f"Re_x {rec[:,0].min():.3g}..{rec[:,0].max():.3g})")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    fig = sys.argv[1]
    ex = sys.argv[2] if len(sys.argv) > 2 else "../exptData/T3A.dat"
    main(fig, ex)
