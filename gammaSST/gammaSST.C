/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | gammaSST-of13 community contribution
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of a community extension for OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gammaSST.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
gammaSST<BasicMomentumTransportModel>::TuL() const
{
    // Local turbulence intensity [%], Menter et al. (2015) Eq. (9):
    //     Tu_L = min(100*sqrt(2 k/3)/(omega y), 100)
    const volScalarField::Internal& k = this->k_();
    const volScalarField::Internal& omega = this->omega_();
    const volScalarField::Internal& y = this->y()();

    // Stabilised against omega*y -> 0 on the wall face itself
    const dimensionedScalar omegaY0
    (
        "small",
        dimVelocity,
        rootVSmall
    );

    return volScalarField::Internal::New
    (
        this->groupName("TuL"),
        min(100*sqrt((2.0/3.0)*k)/max(omega*y, omegaY0), scalar(100))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
gammaSST<BasicMomentumTransportModel>::FPG() const
{
    // Pressure-gradient function, Menter et al. (2015) Eq. (11)-(12).
    //
    //     lambda_thetaL = -7.57e-3 * (dU/dy . n) * y^2/nu + 0.0128
    //     clipped to [-1, 1]
    //
    //     lambda >= 0:  FPG = min(1 + CPG1*lambda, CPG1lim)
    //     lambda <  0:  FPG = min(1 + CPG2*lambda
    //                          + CPG3*min(lambda + 0.0681, 0), CPG2lim)
    //     FPG = max(FPG, 0)
    //
    // n is the wall-normal direction, obtained from the gradient of the
    // wall distance.
    const tmp<volScalarField> tnu = this->nu();
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y()();

    const volVectorField n(fvc::grad(this->y()));

    const volScalarField::Internal dVdy
    (
        (fvc::grad(this->U_ & n)().v() & n.v())
    );

    const volScalarField::Internal lambdaThetaL
    (
        min
        (
            max
            (
                -7.57e-3*dVdy*sqr(y)/nu + 0.0128,
                dimensionedScalar(dimless, -1)
            ),
            dimensionedScalar(dimless, 1)
        )
    );

    tmp<volScalarField::Internal> tFPG
    (
        volScalarField::Internal::New
        (
            this->groupName("FPG"),
            this->mesh_,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField::Internal& FPG = tFPG.ref();

    forAll(FPG, celli)
    {
        const scalar lambda = lambdaThetaL[celli];

        FPG[celli] =
            lambda >= 0
          ?
            min(1 + CPG1_.value()*lambda, CPG1lim_.value())
          :
            min
            (
                1
              + CPG2_.value()*lambda
              + CPG3_.value()*min(lambda + 0.0681, scalar(0)),
                CPG2lim_.value()
            );

        FPG[celli] = max(FPG[celli], scalar(0));
    }

    return tFPG;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
gammaSST<BasicMomentumTransportModel>::ReThetac() const
{
    // Transition onset correlation, Menter et al. (2015) Eq. (8):
    //     ReThetac = CTU1 + CTU2*exp(-CTU3*TuL*FPG)
    return volScalarField::Internal::New
    (
        this->groupName("ReThetac"),
        CTU1_ + CTU2_*exp(-CTU3_*TuL()*FPG())
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> gammaSST<BasicMomentumTransportModel>::Fonset
(
    const volScalarField::Internal& S,
    const volScalarField::Internal& RT
) const
{
    // Transition onset function, Menter et al. (2015) Eq. (4)-(7):
    //     Fonset1 = Re_v/(2.2 ReThetac),   Re_v = y^2 S/nu
    //     Fonset2 = min(Fonset1, 2)
    //     Fonset3 = max(1 - (R_T/3.5)^3, 0)
    //     Fonset  = max(Fonset2 - Fonset3, 0)
    const tmp<volScalarField> tnu = this->nu();
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y()();

    const volScalarField::Internal Rev(sqr(y)*S/nu);

    const volScalarField::Internal Fonset1(Rev/(2.2*ReThetac()));

    const volScalarField::Internal Fonset2(min(Fonset1, scalar(2)));

    const volScalarField::Internal Fonset3(max(1 - pow3(RT/3.5), scalar(0)));

    return volScalarField::Internal::New
    (
        this->groupName("Fonset"),
        max(Fonset2 - Fonset3, scalar(0))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> gammaSST<BasicMomentumTransportModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    // Modified F1 blending function, Menter et al. (2015) Eq. (17):
    //     F1 = max(F1_SST, exp(-(R_y/120)^8)),   R_y = y sqrt(k)/nu
    //
    // This keeps the k-omega branch active through the laminar boundary
    // layer, where the k-epsilon branch would otherwise be selected.
    const volScalarField Ry(this->y()*sqrt(this->k_)/this->nu());

    const volScalarField F3(exp(-pow(Ry/120.0, 8)));

    return max(kOmegaSST<BasicMomentumTransportModel>::F1(CDkOmega), F3);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> gammaSST<BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    // Modified production of k, Menter et al. (2015) Eq. (13) and (15):
    //     Pk -> gamma*Pk + Pk_lim
    //
    // PkLim_ is evaluated in correctGammaInt(), which runs immediately
    // before kOmegaSST::correct() calls this function.
    return
        gammaInt_()*kOmegaSST<BasicMomentumTransportModel>::Pk(G)
      + PkLim_;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
gammaSST<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    // Modified destruction of k, Menter et al. (2015) Eq. (14):
    //     Dk -> min(max(gamma, 0.1), 1)*Dk
    return
        min(max(gammaInt_(), scalar(0.1)), scalar(1))
       *kOmegaSST<BasicMomentumTransportModel>::epsilonByk(F1, F2);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
gammaSST<BasicMomentumTransportModel>::gammaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    kOmegaSST<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity,
        type
    ),

    // Intermittency equation coefficients, Menter et al. (2015) Table 1
    Flength_("Flength", this->coeffDict(), 100),
    ca2_("ca2", this->coeffDict(), 0.06),
    ce2_("ce2", this->coeffDict(), 50),
    sigmaGamma_("sigmaGamma", this->coeffDict(), 1),

    // Transition onset correlation coefficients
    CTU1_("CTU1", this->coeffDict(), 100),
    CTU2_("CTU2", this->coeffDict(), 1000),
    CTU3_("CTU3", this->coeffDict(), 1),

    // Pressure-gradient function coefficients
    CPG1_("CPG1", this->coeffDict(), 14.68),
    CPG1lim_("CPG1lim", this->coeffDict(), 1.5),
    CPG2_("CPG2", this->coeffDict(), -7.34),
    CPG3_("CPG3", this->coeffDict(), 0),
    CPG2lim_("CPG2lim", this->coeffDict(), 3),

    // Separation-induced transition coefficients
    ReThetacLim_("ReThetacLim", this->coeffDict(), 1100),
    Ck_("Ck", this->coeffDict(), 1),
    CSEP_("CSEP", this->coeffDict(), 1),

    gammaInt_
    (
        IOobject
        (
            this->groupName("gammaInt"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    PkLim_
    (
        IOobject
        (
            this->groupName("PkLim"),
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(this->k_.dimensions()/dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool gammaSST<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST<BasicMomentumTransportModel>::read())
    {
        Flength_.readIfPresent(this->coeffDict());
        ca2_.readIfPresent(this->coeffDict());
        ce2_.readIfPresent(this->coeffDict());
        sigmaGamma_.readIfPresent(this->coeffDict());

        CTU1_.readIfPresent(this->coeffDict());
        CTU2_.readIfPresent(this->coeffDict());
        CTU3_.readIfPresent(this->coeffDict());

        CPG1_.readIfPresent(this->coeffDict());
        CPG1lim_.readIfPresent(this->coeffDict());
        CPG2_.readIfPresent(this->coeffDict());
        CPG3_.readIfPresent(this->coeffDict());
        CPG2lim_.readIfPresent(this->coeffDict());

        ReThetacLim_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
        CSEP_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void gammaSST<BasicMomentumTransportModel>::correctGammaInt()
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const tmp<volScalarField> tnu = this->nu();
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y()();

    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    // Strain rate and vorticity magnitudes
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal S(sqrt(2*magSqr(symm(tgradU()()))));
    const volScalarField::Internal Omega(sqrt(2*magSqr(skew(tgradU()()))));
    tgradU.clear();

    // Viscosity ratio, R_T = k/(nu omega)
    const volScalarField::Internal RT(k()/(nu*omega()));

    {
        // Production of intermittency, Menter et al. (2015) Eq. (2):
        //     P_gamma = Flength rho S gamma (1 - gamma) F_onset
        //
        // The leading gamma is folded into Pgamma so that the bracket left
        // for the implicit/explicit split below is the (1 - gamma) factor.
        const volScalarField::Internal Pgamma
        (
            alpha()*rho()*Flength_*S*Fonset(S, RT)*gammaInt_()
        );

        // Destruction of intermittency, Menter et al. (2015) Eq. (3):
        //     E_gamma = ca2 rho Omega gamma F_turb (ce2 gamma - 1)
        //
        // Likewise the leading gamma is folded into Egamma.
        const volScalarField::Internal Fturb(exp(-pow4(RT/2)));

        const volScalarField::Internal Egamma
        (
            alpha()*rho()*ca2_*Omega*Fturb*gammaInt_()
        );

        // Intermittency equation.
        //
        // Source terms use the standard OpenFOAM implicit/explicit split:
        // the gamma-linear part of each term goes into the matrix diagonal
        // via fvm::Sp.  Note that BOTH Sp contributions are NEGATIVE on the
        // right-hand side, which is what keeps the matrix diagonally
        // dominant; a positive Sp here makes the equation unstable and lets
        // gamma run away above 1.
        //
        //   production   Pgamma - Sp(Pgamma, gamma)     ->  Pgamma (1 - gamma)
        //   destruction  Egamma - Sp(ce2 Egamma, gamma) ->  Egamma (1 - ce2 gamma)
        //                                               =  -E_gamma  (as required,
        //                                                  since the equation is
        //                                                  d(rho gamma)/dt = P - E)
        //
        // This matches the sign convention used by the built-in kOmegaSSTLM.
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
          + fvm::div(alphaRhoPhi, gammaInt_)
          - fvm::laplacian(alpha*rho*DgammaIntEff(), gammaInt_)
         ==
            Pgamma - fvm::Sp(Pgamma, gammaInt_)
          + Egamma - fvm::Sp(ce2_*Egamma, gammaInt_)
          + fvModels.source(alpha, rho, gammaInt_)
        );

        gammaIntEqn.ref().relax();
        fvConstraints.constrain(gammaIntEqn.ref());
        solve(gammaIntEqn);
        fvConstraints.constrain(gammaInt_);
        bound(gammaInt_, 0);
    }

    // Separation-induced transition correction, Menter et al. (2015)
    // Eq. (15)-(16):
    //     Fon_lim = min(max(Re_v/(2.2 ReThetac_lim) - 1, 0), 3)
    //     Pk_lim  = 5 Ck max(gamma - 0.2, 0) (1 - gamma) Fon_lim
    //               max(3 CSEP nu - nut, 0) S Omega
    const volScalarField::Internal Rev(sqr(y)*S/nu);

    const volScalarField::Internal FonLim
    (
        min(max(Rev/(2.2*ReThetacLim_) - 1, scalar(0)), scalar(3))
    );

    PkLim_ =
        5*Ck_
       *max(gammaInt_() - 0.2, scalar(0))
       *(1 - gammaInt_())
       *FonLim
       *max(3*CSEP_*nu - this->nut_(), dimensionedScalar(nu.dimensions(), 0))
       *S*Omega;
}


template<class BasicMomentumTransportModel>
void gammaSST<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Solve the intermittency equation and evaluate the separation-induced
    // correction to the production of k
    correctGammaInt();

    // Solve k and omega; the overridden F1(), Pk() and epsilonByk() apply
    // the transition modifications
    kOmegaSST<BasicMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
