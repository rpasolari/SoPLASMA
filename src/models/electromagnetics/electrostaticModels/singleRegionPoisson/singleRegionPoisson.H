/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::singleRegionPoisson
    
Description
    Poisson equation solver for the electrostatic potential in a single
    mesh region (gas only, no dielectrics).

    Solves:
        ∇·(ε∇φ) = -ρ_e                          (explicit)
        ∇·((ε + Δt·σ)∇φ) = -ρ_e - Δt·S_diff    (semi-implicit)

    chargeDensity and surfCharge are inherited from the base class. The
    calling solver updates them via mutable references before each solve.
    Conductivity and diffusive source are passed as arguments since they
    are purely transport quantities with no reason to persist.

SourceFiles
    singleRegionPoisson.C

\*---------------------------------------------------------------------------*/

#ifndef singleRegionPoisson_H
#define singleRegionPoisson_H

#include "electromagneticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class singleRegionPoisson Declaration
\*---------------------------------------------------------------------------*/

class singleRegionPoisson
:
    public electromagneticsModel
{
private:

    //- Disallow copy
    singleRegionPoisson(const singleRegionPoisson&) = delete;
    void operator=(const singleRegionPoisson&) = delete;

    //- Update E, Emag, phiE, and reducedE from ePotential
    void updateDerivedFields();

protected:

    //- Numerical scheme for computing E
    word EScheme_;

    //- Numerical scheme (explicit vs semiImplicit)
    word PoissonScheme_;

    //- Number of non-orthogonal corrector passes
    label nNonOrthCorr_;

    //- Uniform background gas density for reducedE calculations [1/m3]
    dimensionedScalar backgroundDensityUniform_;

    //- Pointer to spatially varying background density
    const volScalarField* backgroundDensityFieldPtr_;

public:

    //- Runtime type information
    TypeName("singleRegionPoisson");

    // Constructors

        //- Construct from meshes
        singleRegionPoisson
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    //- Destructor
    virtual ~singleRegionPoisson() = default;

    // Solve Functions

        //- Solve Poisson equation (for explicit Poisson)
        virtual void solve() override;

        //- Solve semi-implicit Poisson equation (for semiImplicit Poisson)
        virtual void solve
        (
            const volScalarField& electricalConductivity,
            const volScalarField& diffusiveChargeSource
        ) override;

    // Configuration

        //- Set uniform background density for reducedE
        void setBackgroundDensity(const dimensionedScalar& N)
        {
            backgroundDensityUniform_ = N;
            backgroundDensityFieldPtr_ = nullptr;
        }

        //- Set spatially varying background density for reducedE
        void setBackgroundDensity(const volScalarField& N)
        {
            backgroundDensityFieldPtr_ = &N;
        }

        //- Return the Poisson scheme
        const word& PoissonScheme() const override { return PoissonScheme_; }

        //- Return number of non-orthogonal corrector passes
        label nNonOrthCorr() const override { return nNonOrthCorr_; }

        //- Return E-field scheme
        const word& EScheme() const override { return EScheme_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

