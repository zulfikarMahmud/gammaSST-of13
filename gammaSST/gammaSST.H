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

Class
    Foam::RASModels::gammaSST

Description
    Menter-Smirnov one-equation local correlation-based transition model
    (the "gamma" model), coupled to the k-omega-SST turbulence model.

    A single intermittency (gamma) transport equation replaces the two
    transport equations of the Langtry-Menter gamma-ReThetat model: the
    transition onset momentum-thickness Reynolds number is instead obtained
    from an algebraic correlation of the local turbulence intensity and
    pressure-gradient parameter, which also makes the model Galilean
    invariant (no dependence on a free-stream reference velocity).

    Reference:
    \verbatim
        Menter, F. R., Smirnov, P. E., Liu, T., & Avancha, R. (2015).
        A one-equation local correlation-based transition model.
        Flow, Turbulence and Combustion, 95(4), 583-619.
        doi:10.1007/s10494-015-9622-4
    \endverbatim

    The model is implemented exactly as published: no additional terms and
    no modified coefficients.  The governing equation is

    \f[
        \frac{\partial (\rho \gamma)}{\partial t}
      + \frac{\partial (\rho U_j \gamma)}{\partial x_j}
      = P_\gamma - E_\gamma
      + \frac{\partial}{\partial x_j}
        \left[
            \left( \mu + \frac{\mu_t}{\sigma_\gamma} \right)
            \frac{\partial \gamma}{\partial x_j}
        \right]
    \f]

    and it modifies the underlying k-omega-SST model in three places:
    the production and destruction of k, and the F1 blending function.

    The default model coefficients are those of the original paper
    \verbatim
        gammaSSTCoeffs
        {
            // Default SST coefficients are inherited from kOmegaSST

            // Intermittency equation (Menter et al. 2015)
            Flength     100;
            ca2         0.06;
            ce2         50;
            sigmaGamma  1;

            // Transition onset correlation, ReThetac(TuL, FPG)
            CTU1        100;
            CTU2        1000;
            CTU3        1;

            // Pressure-gradient function FPG(lambdaThetaL)
            CPG1        14.68;
            CPG1lim     1.5;
            CPG2        -7.34;
            CPG3        0;
            CPG2lim     3;

            // Separation-induced transition correction
            ReThetacLim 1100;
            Ck          1;
            CSEP        1;
        }
    \endverbatim

    Boundary conditions for the intermittency field \c gammaInt are
    \c inletOutlet (or \c fixedValue) with a value of 1 at inlets and
    \c zeroGradient at walls.

SourceFiles
    gammaSST.C

\*---------------------------------------------------------------------------*/

#ifndef gammaSST_H
#define gammaSST_H

// The kOmegaSST base-class template implementation (kOmegaSSTBase.C, pulled
// in by kOmegaSSTBase.H under NoRepository) uses the fvm:: operators but does
// not include their headers itself.  In the core library they happen to be
// available because an earlier model in the same translation unit includes
// eddyViscosity.H.  Include them explicitly here so that this header is
// self-contained and can be compiled on its own.
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcMeshPhi.H"

#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class gammaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class gammaSST
:
    public kOmegaSST<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Model coefficients

            //- Transition length scale constant
            dimensionedScalar Flength_;

            //- Intermittency destruction coefficient
            dimensionedScalar ca2_;

            //- Intermittency destruction coefficient
            dimensionedScalar ce2_;

            //- Turbulent Prandtl number for the intermittency equation
            dimensionedScalar sigmaGamma_;

            //- Coefficients of the ReThetac(TuL, FPG) correlation
            dimensionedScalar CTU1_;
            dimensionedScalar CTU2_;
            dimensionedScalar CTU3_;

            //- Coefficients of the pressure-gradient function FPG
            dimensionedScalar CPG1_;
            dimensionedScalar CPG1lim_;
            dimensionedScalar CPG2_;
            dimensionedScalar CPG3_;
            dimensionedScalar CPG2lim_;

            //- Critical Reynolds number limit for the separation correction
            dimensionedScalar ReThetacLim_;

            //- Separation-induced transition coefficients
            dimensionedScalar Ck_;
            dimensionedScalar CSEP_;


        // Fields

            //- Intermittency
            volScalarField gammaInt_;

            //- Separation-induced correction to the production of k
            //  Stored between correct() and the Pk() call made by the
            //  underlying kOmegaSST::correct()
            volScalarField::Internal PkLim_;


    // Protected Member Functions

        //- Local turbulence intensity [%], Menter et al. (2015) Eq. (9)
        tmp<volScalarField::Internal> TuL() const;

        //- Pressure-gradient function, Menter et al. (2015) Eq. (11)
        tmp<volScalarField::Internal> FPG() const;

        //- Critical momentum-thickness Reynolds number correlation,
        //  Menter et al. (2015) Eq. (8)
        tmp<volScalarField::Internal> ReThetac() const;

        //- Transition onset function, Menter et al. (2015) Eq. (4)-(7)
        tmp<volScalarField::Internal> Fonset
        (
            const volScalarField::Internal& S,
            const volScalarField::Internal& RT
        ) const;

        //- Modified form of the k-omega SST F1 blending function,
        //  Menter et al. (2015) Eq. (17)
        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

        //- Modified k production rate: gamma*Pk + PkLim,
        //  Menter et al. (2015) Eq. (13) and (15)
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G
        ) const;

        //- Modified k destruction rate: max(gamma, 0.1)*epsilonByk,
        //  Menter et al. (2015) Eq. (14)
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField::Internal& F1,
            const volScalarField::Internal& F2
        ) const;

        //- Solve the intermittency equation and update PkLim
        void correctGammaInt();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("gammaSST");


    // Constructors

        //- Construct from components
        gammaSST
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        gammaSST(const gammaSST&) = delete;


    //- Destructor
    virtual ~gammaSST()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to the intermittency
        const volScalarField& gammaInt() const
        {
            return gammaInt_;
        }

        //- Return the effective diffusivity for the intermittency
        tmp<volScalarField> DgammaIntEff() const
        {
            return volScalarField::New
            (
                "DgammaIntEff",
                this->nut_/sigmaGamma_ + this->nu()
            );
        }

        //- Solve the transport equations and correct the eddy viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gammaSST&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gammaSST.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
