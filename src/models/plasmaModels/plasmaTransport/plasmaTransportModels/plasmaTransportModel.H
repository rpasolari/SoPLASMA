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
    Foam::plasmaTransportModel

Description
    Abstract base class for per-species plasma transport behaviour.

    This class defines the minimal interface required for describing how a
    single plasma species is transported within the simulation.

    The plasmaTransportModel does not implement any transport physics itself.
    Instead, it serves as a common framework from which concrete transport
    models are derived. Each derived class implements its own formulation
    for updating the transport state of a species during the solution loop.

    The run-time selection mechanism is used to construct the appropriate
    model for each species based on the 'transportModel' entry within the
    species dictionary.

SourceFiles
    plasmaTransportModelI.H
    plasmaTransportModel.C
    
\*---------------------------------------------------------------------------*/

#ifndef plasmaTransportModel_H
#define plasmaTransportModel_H

#include "fvMesh.H"
#include "fvMatrix.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

#include "plasmaSpecies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class plasmaTransportModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaTransportModel
{
private:

    //- Disallow default bitwise copy construct
    plasmaTransportModel(const plasmaTransportModel&) = delete;
    void operator=(const plasmaTransportModel&) = delete;

protected:

    //- Model name (e.g. driftDiffusion, fullMomentum, immobile, etc.)
    const word modelName_;

    //- Mesh reference
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Dictionary for this species
    const dictionary& dict_;

    //- Index of the species
    const label specieIndex_;
    
public:

    //- Runtime type information
    TypeName("plasmaTransportModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,                        // storage class
        plasmaTransportModel,           // base class
        dictionary,                     // lookup keyword
        (                               // constructor signature
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        ),
        (modelName, dict, mesh, species, specieIndex)
    );

    // Constructors

        //- Construct base model
        plasmaTransportModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        );

    //- Selector
    static autoPtr<plasmaTransportModel> New
    (
        const word& modelName,
        const dictionary& dict,
        const fvMesh& mesh,
        plasmaSpecies& species,
        const label specieIndex
    );

    //- Destructor
    virtual ~plasmaTransportModel() = default;

    // Member Functions

        //- Update transport coefficients
        virtual void correct() = 0;

        //- Return the transport equation (Matrix)
        virtual tmp<fvScalarMatrix> nEqn() const = 0;

        //- Update the fluxes in place (1/s)
        virtual void updateFluxes
        (
            const fvScalarMatrix& nEqnMatrix,
            surfaceScalarField& convectiveFlux,
            surfaceScalarField& diffusiveFlux,
            surfaceScalarField& particleFlux
        ) const = 0;

        //- Return electrical conductivity contribution [S/m]
        virtual tmp<volScalarField> electricalConductivity() const = 0;

        //- Return diffusive charge source contribution [C/m3/s]
        virtual tmp<volScalarField> diffusiveChargeSource() const = 0;

    // Access

        inline const word& modelName() const;
        inline const fvMesh& mesh() const;
        inline const plasmaSpecies& species() const;
        inline label specieIndex() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

