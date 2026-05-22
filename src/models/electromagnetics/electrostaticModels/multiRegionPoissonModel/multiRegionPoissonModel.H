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
    Foam::multiRegionPoissonModel
    
Description
    Poisson equation solver for the electrostatic potential in multiple
    domains (gas + dielectric regions).

    In the gas/plasma region:
        ∇·(ε∇φ) = -ρ_e                          (explicit)
        ∇·((ε + Δt·σ)∇φ) = -ρ_e - Δt·S_diff    (semi-implicit)

    In each dielectric region:
        ∇·(ε∇φ) = 0

    Supports two coupling architectures:
        - Monolithic: all regions in a single global matrix.
          Requires 'useImplicit true' on mapped patches.
        - Segregated: regions solved iteratively until convergence.

    chargeDensity and surfCharge (gas side) are inherited from the base.
    Each dielectricRegion owns its own surfCharge for the dielectric side
    of the interface.

SourceFiles
    multiRegionPoissonModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegionPoissonModel_H
#define multiRegionPoissonModel_H

#include "electromagneticsModel.H"
#include "dielectricRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiRegionPoissonModel Declaration
\*---------------------------------------------------------------------------*/

class multiRegionPoissonModel
:
    public electromagneticsModel
{
private:

    //- Disallow copy
    multiRegionPoissonModel(const multiRegionPoissonModel&) = delete;
    void operator=(const multiRegionPoissonModel&) = delete;

    //- Update E, Emag, phiE, and reducedE from ePotential
    void updateDerivedFields();

    //- Monolithic solve (explicit Poisson))
    void solveCoupled();

    //- Monolithic solve (semi-implicit Poisson)
    void solveCoupled
    (
        const volScalarField& electricalConductivity,
        const volScalarField& diffusiveChargeSource
    );

    //- Segregated solve (explicit Poisson)
    void solveSegregated();

    //- Segregated solve (semi-implicit Poisson)
    void solveSegregated
    (
        const volScalarField& electricalConductivity,
        const volScalarField& diffusiveChargeSource
    );


protected:

    //- Numerical scheme for computing E ("grad" or "reconstruct")
    word EScheme_;

    //- Numerical scheme (explicit vs semiImplicit)
    word PoissonScheme_;

    //- Number of non-orthogonal corrector passes
    label nNonOrthCorr_;

    //- Uniform background gas density for reducedE calculations [1/m3]
    dimensionedScalar backgroundDensityUniform_;

    //- Pointer to spatially varying background density (not owned)
    const volScalarField* backgroundDensityFieldPtr_;

    //- A ptr list for the dielectric regions
    PtrList<dielectricRegion> dielectrics_;

    //- A boolean to check if the main region is coupled with the dielectrics
    Switch coupled_;

    //- Matrix assembly for monolithic solve
    autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr_;

    //- Max iterations for segregated (non-coupled) loop
    label maxNonCoupledIterations_;

    //- Convergence tolerance for (non-coupled) segregated loop
    scalar nonCoupledTolerance_;

public:

    //- Runtime type information
    TypeName("multiRegionPoisson");

    // Constructors

        //- Construct from meshes
        multiRegionPoissonModel
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    //- Destructor
    virtual ~multiRegionPoissonModel() = default;

    // Solve

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

        //- Return true if regions are monolithically coupled
        bool coupled() const { return coupled_; }

        //- Return number of non-orthogonal corrector passes
        label nNonOrthCorr() const override { return nNonOrthCorr_; }

        //- Return E-field scheme
        const word& EScheme() const override { return EScheme_; }

    // Dielectric access
    
        //- Dielectric regions container access
        inline const PtrList<dielectricRegion>& dielectrics() const 
        { 
            return dielectrics_; 
        }
        inline label nDielectrics() const { return dielectrics_.size(); }

        //- Access specific dielectric region (by index or name)
        inline const dielectricRegion& dielectric(const label i) const 
        { 
            return dielectrics_[i]; 
        }
        const dielectricRegion& dielectric(const word& name) const;

        //- Access specific dielectric permittivity (by index or name)
        inline const dimensionedScalar& epsilon(const label i) const 
        { 
            return dielectrics_[i].epsilon(); 
        }
        inline const dimensionedScalar& epsilon(const word& name) const 
        { 
            return dielectric(name).epsilon(); 
        }


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

