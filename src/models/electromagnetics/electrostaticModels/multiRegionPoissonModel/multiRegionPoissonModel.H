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
    A Poisson equation solver for electrostatic potential in a single 
    mesh region.

    Solves the Poisson equation for the electric potential:

        ∇·(ε∇φ) = -ρ_e

    where ε is the permittivity and ρ_e is the charge density.

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
    pure electrostatic Laplace equation (∇²φ = 0) within the region.

SourceFiles
    multiRegionPoissonModel.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegionPoissonModel_H
#define multiRegionPoissonModel_H

#include "electromagneticsModel.H"
#include "plasmaTransportModel.H"

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
    //--- Fields ---
    volScalarField ePotential_;
    volVectorField E_;
    dimensionedScalar epsilon_;

    PtrList<volScalarField> ePotentialDielectricList_;
    PtrList<volVectorField> EDielectricList_;
    PtrList<dimensionedScalar> epsilonDielectricList_;

    //--- Persistent Coupling State (Your Logic) ---
    bool coupled_;
    autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr_;

public:

    TypeName("multiRegionPoissonModel");

    // FULL SIGNATURE
    multiRegionPoissonModel
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    );

    virtual ~multiRegionPoissonModel() {}

    virtual void solve() override;

    virtual const volVectorField& E() const override
    {
        return E_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

