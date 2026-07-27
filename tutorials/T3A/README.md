# ERCOFTAC T3A — flat-plate bypass transition

The standard first validation case for a transition model: a sharp-edged flat
plate in a high-turbulence freestream, where transition is triggered by
freestream turbulence (bypass transition) rather than by Tollmien–Schlichting
waves.

## Conditions

| | |
|---|---|
| Freestream velocity `U∞` | 5.4 m/s |
| Kinematic viscosity `ν` | 1.5e-5 m²/s |
| Inlet turbulence intensity | 3.3 % (`k = 0.047633`, `ω = 264.63`) |
| Pressure gradient | zero |
| Leading edge at | `x = 0.04 m` |
| Mesh | 26 820 cells, `y⁺ < 1` |
| Solver | `foamRun -solver incompressibleFluid` (steady SIMPLE) |

Experimental data: `validation/exptData/T3A.dat` — Savill (1993, 1996),
ERCOFTAC Classic Collection.

## Running

```bash
source /opt/openfoam13/etc/bashrc
./Allrun
```

Then, for numbers rather than a picture:

```bash
cd validation && python3 compare.py
```

and for the figure:

```bash
cd validation && python3 plot.py
```

To overlay a second model (e.g. the built-in 4-equation `kOmegaSSTLM`), copy
OpenFOAM's own T3A tutorial, run it, and point `plot.py`/`compare.py` at it:

```bash
cp -r $FOAM_TUTORIALS/incompressibleFluid/T3A /tmp/T3A-LM
cd /tmp/T3A-LM && blockMesh && foamRun && cd -
cd validation && python3 plot.py --compare /tmp/T3A-LM --label "kOmegaSSTLM (built-in)"
```

## Results

![T3A validation](validation/T3A_validation.png)

Both models were run on the **identical mesh, identical boundary conditions and
identical numerics** — the only difference is the turbulence model — so this is a
clean like-for-like comparison.

| Metric | **gammaSST** (this repo) | `kOmegaSSTLM` (built-in) | Experiment |
|---|---|---|---|
| Transition onset `Re_x` (c_f minimum) | **1.42e5**  (ratio 1.01) | 1.09e5  (ratio 0.77) | 1.41e5 |
| `c_f` mean abs. error | **0.00034** (11.2 %) | 0.00039 (13.7 %) | — |
| Freestream `Tu` decay mean abs. error | 0.096 % | 0.095 % | — |
| Transport equations | **3** | 4 | — |
| SIMPLE iterations to convergence | **251** | 268 | — |

**Reading the result.** The `gammaSST` transition onset lands within **1 %** of
the measured `c_f` minimum, where the built-in 4-equation model fires about 23 %
early. This is the expected and published behaviour — improved onset prediction
under a specified freestream turbulence level is precisely what the algebraic
`Reθc` correlation plus the `F_PG` pressure-gradient function were introduced to
achieve in Menter et al. (2015).

Downstream of transition both models collapse onto the same fully-turbulent
`c_f` curve and track the data closely, and the freestream turbulence decay is
essentially exact for both — confirming the `k`–`ω` side of the model is
untouched and behaving as standard SST.

The intermittency stays inside its physical bounds
(`0.020 ≤ γ ≤ 1.001`; the 0.1 % overshoot is ordinary discretisation, since the
γ-equation self-limits at 1 rather than being clipped there), and `k` peaks at
0.19 m²/s² — the correct order for this case.

> **A note on debugging, kept deliberately.** The first working build of this
> model produced `γ ≈ 6300`, `k ≈ 133`, and no transition at all. The cause was
> the sign of the implicit/explicit split in the γ-equation destruction term:
> writing `- Egamma + fvm::Sp(ce2*Egamma, γ)` adds destruction as *production*
> and puts a positive coefficient on the matrix diagonal. The correct form is
> `+ Egamma - fvm::Sp(ce2*Egamma, γ)`, matching the built-in `kOmegaSSTLM`.
> **`γ` outside `[0, 1]` is the diagnostic that catches this immediately** — which
> is why `compare.py` checks it. See §10.5 of the top-level README.

## Cross-validation against furstj/gammaSST

Checked against [furstj/gammaSST](https://github.com/furstj/gammaSST), the
reference ESI-line implementation, on both code and results.

All 15 coefficient defaults are identical, and `ReThetac`, `Fonset`, `Fturb`,
`TuL`, `FPG` and the γ-equation implicit/explicit split are term-for-term
equivalent. His T3A runs at `U∞ = 5.18 m/s` against this case's `5.4 m/s`, with
identical `Tu = 3.3 %` and `ν_t/ν = 12.0`:

| Feature | This repo (OF13) | furstj (published figure) | Experiment |
|---|---|---|---|
| `c_f` minimum | 0.00228 @ `Re_x` 1.42e5 | ≈0.00215 @ ≈1.5e5 | 0.00210 @ 1.41e5 |
| turbulent peak `c_f` | 0.00465 @ 2.7e5 | ≈0.0046 @ ≈3.0e5 | ≈0.0049 @ 3.0e5 |
| `c_f` at `Re_x` = 5e5 | 0.00423 | ≈0.0042 | 0.0042 |

Laminar level, onset, transition length, turbulent peak and downstream recovery
all agree; the residual offset is consistent with the 4 % difference in `U∞`.

His implementation additionally carries a **crossflow** extension (`FonsetCF`,
`CRSF`) from later work, which is deliberately absent here — this repo implements
the 2015 paper as published.

## Difference from OpenFOAM's own T3A tutorial

This case is a copy of `$FOAM_TUTORIALS/incompressibleFluid/T3A` with exactly
three changes:

1. `constant/momentumTransport` — `model kOmegaSSTLM;` → `model gammaSST;`
2. `system/controlDict` — added `libs ("libgammaSST.so");`
3. `0/ReThetat` — deleted (the 4-equation model's field; unused by gammaSST)

plus the added `validation/compare.py` and `validation/plot.py`.

## References

- MENTER, F. R., SMIRNOV, P. E., LIU, T. & AVANCHA, R. (2015). A one-equation
  local correlation-based transition model. *Flow, Turbulence and Combustion*
  **95**(4), 583–619.
- SAVILL, A. M. (1993). Some recent progress in the turbulence modelling of
  by-pass transition. *Near-wall turbulent flows*, 829–848.
- SAVILL, A. M. (1996). One-point closures applied to transition. In
  *Turbulence and Transition Modelling*, 233–268. Springer.
- [ERCOFTAC Classic Collection, case 020](http://cfd.mace.manchester.ac.uk/ercoftac/)
