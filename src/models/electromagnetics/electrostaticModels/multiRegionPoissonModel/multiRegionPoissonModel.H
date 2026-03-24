/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
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
    A multi-region Poisson equation solver for electrostatic potential in 
    multiple domains (gas and dielectric regions).

    In the gas/plasma region, it solves the Poisson equation:
        ∇·(ε∇φ) = -ρ_e
    In the dielectric regions, it solves the Laplace equation:
        ∇·(ε∇φ) = 0

    where ε is the permittivity and ρ_e is the charge density.

    The model supports two coupling architectures:
        - Monolithic (coupled): Solves all regions simultaneously in a single 
          global matrix. Requires 'useImplicit true' on mapped patches.
        - Segregated (non-coupled): Solves regions iteratively, using 
          residual-based convergence control to synchronize potentials at 
          interfaces.

    The model also supports a semi-implicit numerical scheme. By coupling with 
    a plasmaTransportModel, it includes the electrical conductivity (sigma) in 
    the Laplacian operator:

        ∇·((ε + Δt·σ)∇φ) = -ρ_e - Δt·S_diff

    where:
        - Δt is the time step
        - σ is the electrical conductivity (Σ |q_i| * mu_i * n_i)
        - S_diff is the diffusive charge source (∇·Σ q_i * D_i * ∇n_i)

    This model can be initialized with or without a plasmaSpecies object. 
    If no plasmaSpecies is provided, the model defaults to solving the 
    pure electrostatic Laplace equation (∇²φ = 0) within the gas region.

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

protected:
   
    //- Electric potential field [V]
    volScalarField ePotential_;

    //- Electric field [V/m]
    volVectorField E_;

    //- Magnitude of the electric field [V/m]
    volScalarField Emag_;

    //- Flux of the electric field [V * m2]
    surfaceScalarField phiE_;

    //- Dielectric permittivity [F/m]
    dimensionedScalar epsilon_;

    //- Relative permittivity of the medium
    scalar epsilonR_;

    //- Numerical scheme (explicit vs semiImplicit)
    word PoissonScheme_;

    //- Pointer to transport model (can be nullptr)
    const plasmaTransportModel* transportPtr_;

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

        //- Construct from meshes and species
        multiRegionPoissonModel
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes,
            const plasmaSpecies* species
        );

    //- Destructor
    virtual ~multiRegionPoissonModel() = default;

    // Member Functions

        //- The implementation of the solve function
        virtual void solve() override;

        //- The implementation of the solve function for semiImplicit
        virtual void solve(const plasmaTransportModel& transport) override;

        //- Return the LHS of the Poisson equation
        virtual tmp<fvScalarMatrix> PoissonMatrix() const;

        //- Update the internal transport pointer
        void setTransport(const plasmaTransportModel& transport)
        {
            transportPtr_ = &transport;
        }

    // Access

        const volScalarField& ePotential() const override { return ePotential_;}
        const volVectorField& E() const override { return E_; }
        const volScalarField& Emag() const override { return Emag_; }
        const surfaceScalarField& phiE() const override { return phiE_; }
        const dimensionedScalar& epsilon() const override { return epsilon_; }
        const scalar& epsilonR() const override { return epsilonR_; }

        inline const PtrList<dielectricRegion>& dielectrics() const
        {
            return dielectrics_;
        }
        inline label nDielectrics() const
        {
            return dielectrics_.size();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

