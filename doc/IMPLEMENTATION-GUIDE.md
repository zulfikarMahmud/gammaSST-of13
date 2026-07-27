# Implementing a turbulence model in OpenFOAM 13 — step by step

This is the guide I wish had existed. It walks through building the Menter et al.
(2015) γ transition model from nothing, as an **out-of-tree** library that never
touches the OpenFOAM installation.

The *method* generalises: follow the same eight steps for any RAS model you want
to add to OpenFOAM 13.

**Time:** an afternoon if you're comfortable with C++; a weekend if you're
learning the framework as you go.

---

## Contents

- [Step 0 — Understand what you're building](#step-0--understand-what-youre-building)
- [Step 1 — Choose the architecture: inherit, don't copy](#step-1--choose-the-architecture-inherit-dont-copy)
- [Step 2 — Create the skeleton](#step-2--create-the-skeleton)
- [Step 3 — `Make/files`](#step-3--makefiles)
- [Step 4 — `Make/options`](#step-4--makeoptions)
- [Step 5 — The header](#step-5--the-header)
- [Step 6 — The implementation](#step-6--the-implementation)
- [Step 7 — Register for both regimes](#step-7--register-for-both-regimes)
- [Step 8 — Build](#step-8--build)
- [Step 9 — Set up a case](#step-9--set-up-a-case)
- [Step 10 — Validate (do not skip)](#step-10--validate-do-not-skip)
- [Debugging checklist](#debugging-checklist)
- [Porting to other OpenFOAM versions](#porting-to-other-openfoam-versions)

---

## Step 0 — Understand what you're building

Before any code, write down the equations and map each term to something
OpenFOAM can express. For the γ model:

| Paper | Code |
|---|---|
| `∂(ργ)/∂t + ∂(ρu_jγ)/∂x_j` | `fvm::ddt(alpha, rho, gammaInt_) + fvm::div(alphaRhoPhi, gammaInt_)` |
| `∂/∂x_j[(μ + μ_t/σ_γ)∂γ/∂x_j]` | `fvm::laplacian(alpha*rho*DgammaIntEff(), gammaInt_)` |
| `P_γ = F_length ρ S γ(1−γ) F_onset` | explicit + `fvm::Sp` split (Step 6d) |
| `E_γ = c_a2 ρ Ω γ F_turb (c_e2 γ − 1)` | explicit + `fvm::Sp` split (Step 6d) |
| `S`, `Ω` | `sqrt(2*magSqr(symm(gradU)))`, `sqrt(2*magSqr(skew(gradU)))` |
| wall distance `d` | `this->y()` |
| `R_T = ρk/(μω)` | `k()/(nu*omega())` |

**Also write down which existing model yours is a modification of.** Ours is
k-ω SST with:

- one extra transport equation (`γ`),
- `P_k → γ P_k + P_k^lim`,
- `D_k → min(max(γ,0.1),1) D_k`,
- `F1 → max(F1, F3)`.

That list is the whole design. Everything else is inherited.

---

## Step 1 — Choose the architecture: inherit, don't copy

**The single most important decision.** The tempting move is to copy
`kOmegaSST.C` and edit it. Don't. OpenFOAM 13 deliberately exposes the terms you
need as `virtual` functions.

Find them:

```bash
grep -n "virtual" $FOAM_SRC/MomentumTransportModels/momentumTransportModels/Base/kOmegaSST/kOmegaSSTBase.H
```

The relevant hooks:

| Hook | SST default | What we override it to |
|---|---|---|
| `Pk(G)` | `G` | `γ*G + PkLim` |
| `epsilonByk(F1,F2)` | `betaStar*omega` | `min(max(γ,0.1),1) * ...` |
| `F1(CDkOmega)` | SST blending | `max(F1_SST, F3)` |
| `correct()` | solves k, ω | solve γ first, then call base |

**Payoff:** ~480 lines instead of ~1500, and every upstream SST fix propagates
to you for free. If your model modifies a term that is *not* virtual, that is
the one case where you must copy — check first.

Verify the base class name and template signature you must inherit:

```bash
sed -n '1,80p' $FOAM_SRC/MomentumTransportModels/momentumTransportModels/RAS/kOmegaSST/kOmegaSST.H
```

---

## Step 2 — Create the skeleton

Out-of-tree means: anywhere you like, **not** inside `$FOAM_SRC`.

```bash
mkdir -p ~/OpenFOAM/gammaSST-of13/gammaSST/Make
cd ~/OpenFOAM/gammaSST-of13
```

Final layout:

```
gammaSST-of13/
├── Allwmake
├── Allwclean
└── gammaSST/
    ├── Make/{files,options}
    ├── gammaSST.H
    ├── gammaSST.C
    ├── gammaSSTIncompressibleMomentumTransportModels.C
    └── gammaSSTCompressibleMomentumTransportModels.C
```

`Allwmake` is two useful lines:

```sh
#!/bin/sh
cd "${0%/*}" || exit 1
wmake libso gammaSST
```

---

## Step 3 — `Make/files`

Tells `wmake` what to compile and **where the library goes**:

```make
gammaSSTIncompressibleMomentumTransportModels.C
gammaSSTCompressibleMomentumTransportModels.C

LIB = $(FOAM_USER_LIBBIN)/libgammaSST
```

Two things people get wrong here:

1. **`$(FOAM_USER_LIBBIN)`, not `$(FOAM_LIBBIN)`.** The former is
   `~/OpenFOAM/<user>-13/platforms/.../lib` — your own directory. The latter
   writes into the OpenFOAM installation, which is what we're avoiding.
2. **`gammaSST.C` is deliberately absent.** It's a template; it gets `#include`d
   by the two registration units. Compiling it standalone produces no symbols.

---

## Step 4 — `Make/options`

Don't guess the include paths — copy them from the framework you're extending:

```bash
cat $FOAM_SRC/MomentumTransportModels/incompressible/Make/options
cat $FOAM_SRC/MomentumTransportModels/compressible/Make/options
```

Take the union, plus the shared base:

```make
EXE_INC = \
    -I$(LIB_SRC)/MomentumTransportModels/momentumTransportModels/lnInclude \
    -I$(LIB_SRC)/MomentumTransportModels/incompressible/lnInclude \
    -I$(LIB_SRC)/MomentumTransportModels/compressible/lnInclude \
    -I$(LIB_SRC)/physicalProperties/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude

LIB_LIBS = \
    -lmomentumTransportModels \
    -lincompressibleMomentumTransportModels \
    -lcompressibleMomentumTransportModels \
    -lphysicalProperties \
    -lfiniteVolume \
    -lmeshTools
```

`-lphysicalProperties` is needed by the compressible side. Leaving it out
compiles fine and fails at **link** time with undefined symbols.

---

## Step 5 — The header

```cpp
template<class BasicMomentumTransportModel>
class gammaSST
:
    public kOmegaSST<BasicMomentumTransportModel>
{
protected:
    // Model coefficients -- every one from the paper, all runtime-readable
    dimensionedScalar Flength_;
    dimensionedScalar ca2_;
    dimensionedScalar ce2_;
    dimensionedScalar sigmaGamma_;
    dimensionedScalar CTU1_, CTU2_, CTU3_;
    dimensionedScalar CPG1_, CPG1lim_, CPG2_, CPG3_, CPG2lim_;
    dimensionedScalar ReThetacLim_, Ck_, CSEP_;

    // Fields
    volScalarField gammaInt_;              // the transported intermittency
    volScalarField::Internal PkLim_;       // separation-induced correction

    // Correlations
    tmp<volScalarField::Internal> TuL() const;
    tmp<volScalarField::Internal> FPG() const;
    tmp<volScalarField::Internal> ReThetac() const;
    tmp<volScalarField::Internal> Fonset(...) const;

    // The inherited hooks we override
    virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
    virtual tmp<volScalarField::Internal> Pk(const volScalarField::Internal& G) const;
    virtual tmp<volScalarField::Internal> epsilonByk(...) const;

    void correctGammaInt();

public:
    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField   rhoField;

    TypeName("gammaSST");        // <-- the name used in momentumTransport

    // ... constructor ...

    virtual bool read();
    virtual void correct();
};
```

### Two OF13 idioms to internalise

**`volScalarField::Internal`** — cell-centre values only, no boundary field.
Source terms use it: it's cheaper, and it's what the `Pk`/`epsilonByk`
signatures require. Convert from a full field with `field()` or `field.v()`.

**`tmp<...>`** — reference-counted temporary. Returning `tmp<volScalarField>`
avoids copying an entire field. Call `.ref()` for a writable reference,
`()` to dereference.

**`TypeName("gammaSST")`** is what makes `model gammaSST;` resolve at run time.

---

## Step 6 — The implementation

Write `gammaSST.C` in dependency order.

### 6a. Correlations

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
gammaSST<BasicMomentumTransportModel>::ReThetac() const
{
    return CTU1_ + CTU2_*exp(-CTU3_*TuL()*FPG());
}
```

This one algebraic line is what replaces the entire fourth transport equation of
the older γ–Reθt model.

### 6b. The overridden hooks

```cpp
// Pk -> gamma*Pk + Pk_lim
return gammaInt_()*kOmegaSST<BasicMomentumTransportModel>::Pk(G) + PkLim_;

// Dk -> min(max(gamma, 0.1), 1)*Dk
return min(max(gammaInt_(), scalar(0.1)), scalar(1))
      *kOmegaSST<BasicMomentumTransportModel>::epsilonByk(F1, F2);

// F1 -> max(F1_SST, F3)
const volScalarField Ry(this->y()*sqrt(this->k_)/this->nu());
const volScalarField F3(exp(-pow(Ry/120.0, 8)));
return max(kOmegaSST<BasicMomentumTransportModel>::F1(CDkOmega), F3);
```

Note each one **calls the base implementation** and modifies its result. That is
the whole point of Step 1.

### 6c. `read()`

Lets users retune coefficients without recompiling:

```cpp
if (kOmegaSST<BasicMomentumTransportModel>::read())
{
    Flength_.readIfPresent(this->coeffDict());
    // ... one line per coefficient ...
    return true;
}
return false;
```

### 6d. `correct()` — and the source-term split

Order matters: **γ must be solved before the base solves k**, because `Pk()`
reads `gammaInt_`.

```cpp
void gammaSST<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_) return;
    correctGammaInt();                                  // solves gamma, sets PkLim_
    kOmegaSST<BasicMomentumTransportModel>::correct();  // solves k and omega
}
```

The γ equation itself:

```cpp
// Both source terms carry their leading gamma factor
const volScalarField::Internal Pgamma
(
    alpha()*rho()*Flength_*S*Fonset(S, RT)*gammaInt_()
);
const volScalarField::Internal Egamma
(
    alpha()*rho()*ca2_*Omega*Fturb*gammaInt_()
);

tmp<fvScalarMatrix> gammaIntEqn
(
    fvm::ddt(alpha, rho, gammaInt_)
  + fvm::div(alphaRhoPhi, gammaInt_)
  - fvm::laplacian(alpha*rho*DgammaIntEff(), gammaInt_)
 ==
    Pgamma - fvm::Sp(Pgamma, gammaInt_)
  + Egamma - fvm::Sp(ce2_*Egamma, gammaInt_)
);
```

> ### ⚠ The bug that will bite you
>
> **This is where I lost the most time, so read it carefully.**
>
> `fvm::Sp(a, ψ)` on the right-hand side means `+a·ψ` treated implicitly. So:
>
> - `Pgamma - fvm::Sp(Pgamma, γ)` → `Pgamma·(1 − γ)` ✔
> - `+ Egamma - fvm::Sp(ce2·Egamma, γ)` → `Egamma·(1 − ce2·γ)` = **−E_γ** ✔
>
> I originally wrote `- Egamma + fvm::Sp(ce2*Egamma, γ)`, which is `+E_γ` — it
> **adds destruction as production** and puts a *positive* coefficient on the
> matrix diagonal, destroying diagonal dominance.
>
> The result compiled, ran, and "converged". But γ climbed to **6312**,
> `k` to 133 m²/s², and there was no transition at all.
>
> **Two rules that catch this class of bug:**
> 1. On the RHS, source `Sp` terms should be **negative** — `- fvm::Sp(positive, ψ)`.
>    A positive `Sp` is almost always wrong.
> 2. **Check your bounded quantities.** γ is an intermittency: it *must* lie in
>    [0, 1]. Assert it. Anything outside means the signs or the split are wrong.
>
> Cross-check your split against a model already in the tree that has the same
> structure — `kOmegaSSTLM.C` was the reference that settled it.

### 6e. Numerical safety

Three things must be guarded or you'll get floating-point exceptions:

```cpp
// 1. division -- floor the denominator
k()/max(nu*omega(), dimensionedScalar(dims, SMALL))

// 2. exp() of a large positive -- clip the argument
// 3. pow() of a negative base -- clamp with max(..., 0)
```

Cold starts are the danger: uniform fields, `nut = 0`, possibly `U = 0`
somewhere. The built-in `kOmegaSSTLM` has an unguarded `1/sqr(mag(U))` that
SIGFPEs for exactly this reason — don't repeat it.

**Dimensions are checked at run time.** If a source term's dimensions don't
match the equation, OpenFOAM aborts with a clear message. That catches most
algebra slips for free — lean on it.

---

## Step 7 — Register for both regimes

Pure boilerplate, and it's what gives automatic compressible/incompressible
selection.

`gammaSSTIncompressibleMomentumTransportModels.C`:

```cpp
#include "IncompressibleMomentumTransportModel.H"
#include "incompressibleMomentumTransportModel.H"
#include "makeIncompressibleMomentumTransportModel.H"
#include "gammaSST.H"
#include "gammaSST.C"        // template definition must be visible here

makeRASModel(gammaSST);
```

`gammaSSTCompressibleMomentumTransportModels.C` is identical with
`Compressible`/`compressible`.

### Why this gives automatic selection

`gammaSST` is a *template* on `BasicMomentumTransportModel`:

- **Incompressible instantiation:** `alpha` and `rho` are `geometricOneField` —
  a compile-time constant `1`. Every `alpha()*rho()*P` collapses to `P` with
  zero runtime cost.
- **Compressible instantiation:** they are real fields.

One source, two instantiations, registered into two separate runtime-selection
tables. The **solver** picks the table matching its regime. The user writes
`model gammaSST;` either way — no flag, no branching.

This is exactly how the built-in models do it:

```bash
grep -n kOmegaSSTLM $FOAM_SRC/MomentumTransportModels/*/[ic]*MomentumTransportModels.C
```

---

## Step 8 — Build

```bash
source /opt/openfoam13/etc/bashrc
./Allwmake
```

Verify **both** regimes registered:

```bash
nm -DC $FOAM_USER_LIBBIN/libgammaSST.so | grep -o "incompressibleMomentumTransportModel" | sort -u
nm -DC $FOAM_USER_LIBBIN/libgammaSST.so | grep -o "compressibleMomentumTransportModel"   | sort -u
```

Both must print. If only one does, a registration unit isn't in `Make/files`.

---

## Step 9 — Set up a case

Three edits, no recompilation:

**`system/controlDict`**
```cpp
libs            ("libgammaSST.so");
```

**`constant/momentumTransport`**
```cpp
simulationType RAS;
RAS
{
    model           gammaSST;
    turbulence      on;
}
```

**`0/gammaInt`** — dimensionless, 1 in the freestream, `zeroGradient` at walls.

The fastest route to a working case is to copy an existing tutorial that uses a
related model and change those three things:

```bash
cp -r $FOAM_TUTORIALS/incompressibleFluid/T3A myCase
# edit the three files above; delete 0/ReThetat (4-equation model's field)
```

---

## Step 10 — Validate (do not skip)

**A model that compiles and converges can still be completely wrong.** Mine was
— see the warning in Step 6d. It converged to a garbage state.

Validate against a case with published experimental data. For a transition
model, T3A is the right first test because transition location is the whole
point.

Check, in order:

1. **Bounded quantities are bounded.** `0 ≤ γ ≤ 1`. This is the cheapest and
   highest-yield check you have.
2. **Magnitudes are physical.** `k` peaking at 133 m²/s² in a 5 m/s flow is
   nonsense; ~0.19 is right.
3. **The qualitative feature exists.** Is there a laminar region, a c_f minimum,
   a transitional rise, a turbulent decay?
4. **Then** compare numbers against experiment.
5. **Then** check grid sensitivity — and report it. Transition onset for this
   model moves ~8 % per mesh doubling, so a single-mesh onset number is not a
   validation result.

Compare against a **sibling model on the same mesh** (here, the built-in
`kOmegaSSTLM`) — that isolates your model from mesh and numerics.

---

## Debugging checklist

| Symptom | Likely cause |
|---|---|
| `Unknown RASModel type gammaSST` | `libs` line missing from `controlDict`, or library not built |
| `cannot open shared object file` | forgot to `source .../etc/bashrc` |
| `Cannot find file "gammaInt"` | missing `0/gammaInt` |
| Undefined symbols at link | a `-l...` missing from `Make/options` |
| A bounded quantity leaves its bounds | source-term **sign** or `fvm::Sp` split is wrong |
| Converges to a physically absurd state | same as above — the solver happily converges to garbage |
| SIGFPE on startup | unguarded division / `pow` / `exp` on a cold start |
| Dimension mismatch abort | good news — read the message, it names the term |

---

## Porting to other OpenFOAM versions

| Target | What changes |
|---|---|
| Foundation ≤ v9 | framework is `TurbulenceModels` not `MomentumTransportModels`; `turbulenceProperties` not `momentumTransport`; `makeTurbulenceModel` macros |
| ESI `vXXXX` (openfoam.com) | as above, plus `BasicTurbulenceModel` template parameter naming and dictionary differences |
| A DES/LES variant | inherit from the DES base instead; the same hooks apply |

The *structure* — Steps 1 through 10 — is identical in all of them. Only the
names of the base classes and macros move.

---

## References

- MENTER, F. R., SMIRNOV, P. E., LIU, T. & AVANCHA, R. (2015). A one-equation
  local correlation-based transition model. *Flow, Turbulence and Combustion*
  **95**(4), 583–619. [doi:10.1007/s10494-015-9622-4](https://doi.org/10.1007/s10494-015-9622-4)
- OpenFOAM 13 source guide: <https://cpp.openfoam.org/v13/>
- Reference implementation on the ESI line: <https://github.com/furstj/gammaSST>
