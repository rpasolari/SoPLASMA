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
    Foam::electromagneticsModel
    
Description
    Abstract base class for electromagnetic and electrostatic field solvers.

SourceFiles
    electromagneticsModelI.H
    electromagneticsModel.C

\*---------------------------------------------------------------------------*/

#ifndef electromagneticsModel_H
#define electromagneticsModel_H

#include "fvMesh.H"
#include "UPtrList.H"
#include "runTimeSelectionTables.H"

#include "plasmaSpecies.H"
#include "plasmaTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class plasmaSpecies;

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

    //- Pointer to plasma species
    const plasmaSpecies* speciesPtr_;

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
            const UPtrList<fvMesh>& dielectricMeshes,
            const plasmaSpecies* species
        ),
        (mesh, dielectricMeshes, species)
    );


    // Constructors

        //- Construct from meshes and species
        electromagneticsModel
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes,
            const plasmaSpecies* species
        );

        //- Construct from meshes only
        electromagneticsModel
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );

    // Selectors
    
        //- Selector for meshes and species
        static autoPtr<electromagneticsModel> New
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes,
            const plasmaSpecies* species
        );

        //- Selector for meshes only
        static autoPtr<electromagneticsModel> New
        (
            const fvMesh& mesh,
            const UPtrList<fvMesh>& dielectricMeshes
        );


    //- Destructor
    virtual ~electromagneticsModel() = default;


    // Member functions
    
        //- Solve the electromagnetics
        virtual void solve() = 0;

        //- Solve taking into account the transport (for semiImplicit)
        virtual void solve(const plasmaTransportModel& transport) = 0;

        //- Return the LHS of the Poisson equation (for implicit coupling)
        virtual tmp<fvScalarMatrix> PoissonMatrix() const = 0;

    // Virtual getters
        virtual const volScalarField& ePotential() const = 0;
        virtual const volVectorField& E() const = 0;
        virtual const volScalarField& Emag() const = 0;
        virtual const surfaceScalarField& phiE() const = 0;
        virtual const dimensionedScalar& epsilon() const = 0;
        virtual const scalar& epsilonR() const = 0;

    // Access

        //- Access to the main mesh
        inline const fvMesh& mesh() const;

        //- Access to the list of dielectric meshes
        inline const UPtrList<fvMesh>& dielectricMeshes() const;

        //- Access to the dictionary
        inline const dictionary& dict() const;

        //- Check if plasma data is available
        inline bool hasPlasma() const;

        //- Safe access to plasma species (Fatal error if !hasPlasma)
        inline const plasmaSpecies& species() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "electromagneticsModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
