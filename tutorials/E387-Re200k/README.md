# Eppler 387 — laminar separation bubble at Re = 200 000

A low-Reynolds-number airfoil case: at `Re = 200 000` and `α = 0°` the E387
develops a **laminar separation bubble** on the suction side, which reattaches
turbulent. Getting the bubble in the right place is a much harder test of a
transition model than a flat plate, and it is the regime the model is actually
useful for.

Where [T3A](../T3A/) checks bypass transition under high freestream turbulence,
this case checks **separation-induced transition at low freestream turbulence**
(`Tu = 0.1 %`).

## Conditions

| | |
|---|---|
| Airfoil | Eppler 387, chord 1 m, `α = 0°` |
| Freestream velocity `U∞` | 3 m/s |
| Kinematic viscosity `ν` | 1.5e-5 m²/s |
| Reynolds number | 200 000 |
| Freestream turbulence | `Tu = 0.1 %` → `k = 1.35e-5`, `ω = 0.9` |
| Eddy-viscosity ratio at inlet | `ν_t/ν = 1` |
| Mesh | 40 415 cells, `y⁺ < 1` (supplied, gzipped) |
| Solver | `foamRun -solver incompressibleFluid`, steady SIMPLEC |

Inlet turbulence is set by **inlet conditions only** — there are no sustaining
source terms. `ν_t/ν` is anchored to `Tu` (`nutRatio = ratioRef·Tu/TuRef`,
`ratioRef = 10` at `TuRef = 1 %`), which holds the freestream turbulent length
scale fixed. Change `Tu` in `0/include/caseSettings`; everything else is derived.

Experimental data: `validation/exptData/mcghee_e387_cp_alpha0.csv` — McGhee,
Walker & Millard, **NASA TM-4062 (1988)**, Langley Low-Turbulence Pressure
Tunnel, `Re = 200 000`, `α = 0°`.

## Running

```bash
source /opt/openfoam13/etc/bashrc
./Allrun          # 8-way parallel; edit NP in Allrun to match your machine
```

`Allrun` solves 5000 SIMPLE iterations, reconstructs, then calls
`validate_cp.py` to build the mean `Cp` and compare it with the experiment.
To redo just the comparison, with a different averaging window:

```bash
python3 validate_cp.py . 3 2000      # caseDir, Uinf, tStart
```

## Results

![E387 mean Cp](validation/cp_validation.png)

| region | mean \|ΔCp\| | max \|ΔCp\| |
|---|---|---|
| **overall** | **0.0198** | — |
| upper (suction) | 0.0240 | 0.082 |
| lower (pressure) | 0.0156 | 0.068 |

Suction peak: **CFD −0.581 @ `x/c` 0.322** vs **experiment −0.601 @ `x/c` 0.350**.

The bubble pressure plateau (`x/c ≈ 0.6–0.75`) and the reattachment recovery are
both captured.

## Why the mean, and why the uncertainty band

**The run does not converge to a steady state, and that is physical.** The
laminar separation bubble is genuinely unsteady, so a steady SIMPLE solver
limit-cycles instead of converging: here `Cl` swings **~16 % peak-to-peak** over
the last 2000 iterations (`Cd` ~3.5 %). Residuals plateau around 1e-3 and stop
falling.

A single snapshot is therefore not the answer — the **mean over the limit cycle**
is. Forcing deeper "convergence" with first-order schemes would only quench the
real physics.

`validate_cp.py` builds that mean **in post**, by averaging the per-write-time
surface samples in `postProcessing/surfaces/<time>/airfoilSurf.xy` over
`[2000, 5000]`. Doing it in post rather than with a runtime `fieldAverage` keeps
the run to a single phase and lets the window be changed after the fact.

Because the mean comes from a finite number of snapshots (7 here, at
`writeInterval 500`), the plot shows a **±2 SEM band** rather than hiding the
sampling error. Note where it is widest — `x/c` 0.6–0.75, exactly the bubble
region: **the bubble-plateau agreement is the least certain part of the match**
(mean SEM 0.007 `Cp`, max 0.047). For a tighter band, set `writeInterval 100`
to get ~30 snapshots over the same 5000 iterations.

`Cp = pMean / q` with `q = 0.5·U∞² = 4.5` (OpenFOAM `p` is kinematic, `p/ρ`).

## Numerics notes

**Pressure uses GAMG with the `DIC` smoother, deliberately.** The whole
GaussSeidel family (`GaussSeidel`, `symGaussSeidel`, `DICGaussSeidel`) **SIGFPEs
on the first iteration** of a cold start on this mesh, inside
`GAMGSolver::scale()`: from a uniform `p = 0` field they produce an exactly-zero
coarse-level correction, and `scale()` divides by its norm. `DIC` does not, and
converges the pressure in 7–9 iterations. (PCG/DIC also works but needs ~130.)

**Bounding.** The run logs a few thousand `bounding gammaInt` and `bounding k`
messages. These are `linearUpwind` under/overshoots being clipped, not
divergence: `γ` stays within `[−0.015, 1.003]` and `k` peaks at 0.23 m²/s²,
both physical. If you would rather not carry the risk, switching
`div(phi,gammaInt)` and `div(phi,k)` to `limitedLinear 1` quiets it.

## The mesh

`constant/polyMesh` is **supplied**, not generated — there is no `blockMesh`
step. It is stored gzipped (3 MB instead of 9 MB); OpenFOAM reads `.gz` mesh
files natively, so no manual decompression is needed.

## Files

| File | Purpose |
|---|---|
| `0/include/caseSettings` | the only file to edit — `Tu`, `Uinf`, `ν` and all derived inlet values |
| `0/gammaInt` | intermittency, 1 in freestream, `zeroGradient` at the wall |
| `constant/momentumTransport` | selects `gammaSST`, no coefficient overrides |
| `system/controlDict` | loads `libgammaSST.so`; `surfaces` writes the samples that are averaged |
| `validate_cp.py` | forms `pMean` in post, compares with McGhee, writes the figure |

## References

- McGHEE, R. J., WALKER, B. S. & MILLARD, B. F. (1988). *Experimental Results
  for the Eppler 387 Airfoil at Low Reynolds Numbers in the Langley
  Low-Turbulence Pressure Tunnel.* NASA TM-4062.
- MENTER, F. R., SMIRNOV, P. E., LIU, T. & AVANCHA, R. (2015). A one-equation
  local correlation-based transition model. *Flow, Turbulence and Combustion*
  **95**(4), 583–619.
- MENTER, F. R., LANGTRY, R. & VÖLKER, S. (2006). Transition modelling for
  general purpose CFD codes. *Flow, Turbulence and Combustion* **77**, 277–303.
  (inlet `μ_t/μ` guidance used for the `nutRatio(Tu)` anchoring)
