# Reference data — furstj/gammaSST

## Why this directory exists

[furstj/gammaSST](https://github.com/furstj/gammaSST) is the de-facto reference
implementation of the Menter et al. (2015) γ model, and the obvious thing to
check this port against. Two obstacles:

1. His repository publishes **figures only** — no numerical CFD output is
   committed for the T3A case.
2. His code targets the **ESI** OpenFOAM line (`TurbulenceModels`,
   `turbulentTransportModels`). It cannot be compiled against OpenFOAM 13's
   `MomentumTransportModels` framework without porting it — and porting it is
   precisely what this repository *is*, so building it to "check" this port
   would be circular anyway.

So his `c_f` curve is recovered from his published figure by
`digitise_furstj.py`, giving numbers that can be compared point-by-point.

## The figure is not redistributed here

`digitise_furstj.py` takes a path to his PNG. Clone his repo to regenerate:

```bash
git clone https://github.com/furstj/gammaSST.git
python3 digitise_furstj.py gammaSST/testCases/T3A/validation/figures/Rex_vs_cf.png
```

## Why the digitisation can be trusted

It is **self-calibrating**. His figure plots the experimental points as well as
his curve, and we hold that experimental data exactly — his
`exptData/T3A.dat` is **byte-identical** to OpenFOAM 13's own copy. So the
script predicts where each experimental marker must land from the axis mapping,
detects the actual markers, and compares:

```
calibration on 16 known experimental markers:
    mean 0.44 px, max 0.82 px
    1 px = 761 Re_x, 1.38e-05 c_f
```

Sub-pixel agreement on 16 independent points confirms the axis mapping before
any curve value is read. The script aborts if calibration exceeds 2 px.

Digitisation uncertainty is therefore **±~800 in `Re_x`** and **±1.4e-05 in
`c_f`** — small compared with the differences discussed below.

## Files

| File | What it is |
|---|---|
| `digitise_furstj.py` | the extractor + calibration check |
| `furstj_T3A_cf.csv` | recovered `Re_x, c_f` (734 points), header records provenance |

## Comparison conditions

His T3A runs at `U∞ = 5.18 m/s`; OpenFOAM 13's tutorial uses `5.4 m/s`. Both use
`Tu = 3.3 %` and `ν_t/ν = 12`, and — verified by diff — the **coarse mesh is
identical vertex-for-vertex** between the two repositories.

To remove `U∞` as a variable, the four-way comparison re-runs this model *and*
the built-in `kOmegaSSTLM` at **`U∞ = 5.18 m/s`**, matching him exactly:

```bash
# from tutorials/T3A/validation
python3 plot_compare.py --ours <ourCase> --lm <lmCase> --uinf 5.18
```

See [`../../README.md`](../../README.md) for the results table and discussion.
