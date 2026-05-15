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
    Foam::electromagneticsModel
    
Description
    Abstract base class for electromagnetic and electrostatic field solvers.

    Owns all fields common to plasma EM models: electric potential, electric
    field, charge density, and surface charge. These are persistent, registered
    in the mesh registry, and read from disk where present.

    For pure electrostatic runs, chargeDensity and surfCharge are zero or
    read from disk and never modified. For plasma runs, the calling solver
    updates them via mutable references before each solve.

    The class has no knowledge of plasma species or transport models.
    Conductivity and diffusive sources are passed as raw fields by the solver.

SourceFiles
    electromagneticsModelI.H
    electromagneticsModel.C

\*---------------------------------------------------------------------------*/

#ifndef electromagneticsModel_H
#define electromagneticsModel_H

#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvc.H"
#include "fvm.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "UPtrList.H"
#include "runTimeSelectionTables.H"

#include "SoPLASMAConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class electromagneticsModel Declaration
\*---------------------------------------------------------------------------*/

class electromagneticsModel
:
    public IOdictionary
{
private:

    //- Disallow copy
    electromagneticsModel(const electromagneticsModel&) = delete;
    void operator=(const electromagneticsModel&) = delete;

protected:

    //- Reference to the primary gas/plasma mesh
    const fvMesh& mesh_;

    //- Reference to the list of dielectric meshes
    const UPtrList<fvMesh>& dielectricMeshes_;

    //- Electric potential field [V]
    volScalarField ePotential_;

    //- Electric field [V/m]
    volVectorField E_;

    //- Magnitude of the electric field [V/m]
    volScalarField Emag_;

    //- Flux of the electric field [V * m2]
    surfaceScalarField phiE_;

    //- Reduced electric field [V * m2]
    volScalarField reducedE_;

    //- Voumetric charge density in the gas region [C/m3]
    volScalarField chargeDensity_;

    //- Surface charge density (owned by gas) [C/m2]
    volScalarField surfCharge_;

    //- Dielectric permittivity in the gas region [F/m]
    dimensionedScalar epsilon_;

    //- Relative permittivity of the medium (gas region)
    scalar epsilonR_;

public:

    //- Runtime type information
    TypeName("electromagneticsModel");

    //- Declare the Run-Time Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        electromagneticsModel,
        dictionary,
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        ),
        (mesh, dielectricMeshes)
    );


    // Constructors

        //- Construct from meshes
        electromagneticsModel
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    // Selectors
    
        //- Selector
        static autoPtr<electromagneticsModel> New
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    //- Destructor
    virtual ~electromagneticsModel() = default;

    // Solve functions
    
        //- Solve Poisson equation (for explicit Poisson)
        virtual void solve() = 0;

        //- Solve semi-implicit Poisson equation (for semiImplicit Poisson)
        virtual void solve
        (
            const volScalarField& electricalConductivity,
            const volScalarField& diffusiveChargeSource
        ) = 0;

    // Matrix access
    
        //- Return the LHS Poisson matrix: fvm::laplacian(ε, φ) (explicit)
        virtual tmp<fvScalarMatrix> PoissonLHSMatrix() const = 0;

        //- Return the Poisson equation (explicit Poisson)
        virtual tmp<fvScalarMatrix> PoissonEquation() const = 0;

        //- Return the semi-implicit Poisson equation (semi-implicit Poisson)
        virtual tmp<fvScalarMatrix> PoissonEquation
        (
            const volScalarField& electricalConductivity,
            const volScalarField& diffusiveChargeSource
        ) const = 0;

    // Field access
    
        const volScalarField& ePotential() const { return ePotential_; }
        const volVectorField& E() const { return E_; }
        const volScalarField& Emag() const { return Emag_; }
        const surfaceScalarField& phiE() const { return phiE_; }
        const volScalarField& reducedE() const { return reducedE_; }
        const volScalarField& chargeDensity() const { return chargeDensity_; }
        volScalarField& chargeDensity() { return chargeDensity_; }
        const volScalarField& surfCharge() const { return surfCharge_; }
        volScalarField& surfCharge() { return surfCharge_; }
        const dimensionedScalar& epsilon() const { return epsilon_; }
        scalar epsilonR() const { return epsilonR_; }


    // Mesh and dictionary access

        //- Access to the main mesh
        inline const fvMesh& mesh() const;

        //- Access to the list of dielectric meshes
        inline const UPtrList<fvMesh>& dielectricMeshes() const;

        //- Access to the dictionary
        inline const dictionary& dict() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "electromagneticsModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
