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
    Foam::plasmaEnergyModel

Description
    Abstract base class for per-species plasma energy behaviour.

    This class defines the minimal interface required for describing how the 
    energy in a single plasma species is transported within the simulation.

    The plasmaEnergyModel does not implement any energy physics itself.
    Instead, it serves as a common framework from which concrete energy
    models are derived.  Each derived class implements its own formulation
    for updating the energy state of a species during the solution loop.

    The run-time selection mechanism is used to construct the appropriate
    model for each species based on the 'energyModel' entry within the
    species dictionary.

SourceFiles
    plasmaEnergyModelI.H
    plasmaEnergyModel.C
    
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaEnergyModel_H
#define Foam_plasmaEnergyModel_H

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
                     Class plasmaEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class plasmaEnergyModel
{
protected:

    //- Model name (e.g. isothermal, localField, localEnergy)
    const word modelName_;

    //- Mesh reference
    const fvMesh& mesh_;

    //- Reference to species database
    const plasmaSpecies& species_;

    //- Dictionary for this species
    const dictionary& dict_;

    //- Index of the species this model applies to
    const label specieIndex_;

    //- Reference to the electric field
    const volVectorField& E_;

public:

    //- Runtime type information
    TypeName("plasmaEnergyModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,                        // storage class
        plasmaEnergyModel,           // base class
        dictionary,                     // lookup keyword
        (                               // constructor signature
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        ),
        (modelName, dict, mesh, species, specieIndex, E)
    );

    // Constructors

        //- Construct base model
        plasmaEnergyModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        );

    //- Selector
    static autoPtr<plasmaEnergyModel> New
    (
        const word& modelName,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex,
        const volVectorField& E
    );

    //- Destructor
    virtual ~plasmaEnergyModel() {}

    // Member Functions

        //- Correct model
        virtual void correct() = 0;

        //- Return the energy equation (Matrix)
        virtual tmp<fvScalarMatrix> eEqn() const = 0;

    // Access
    
        //- Returns true if the temperature is uniform (constant in space)
        virtual bool isUniform() const = 0;

        //- Returns true if the temperature is constant in time
        virtual bool isConstant() const = 0;

        //- Return a tmp T field
        //- It can be a reference to a real field or a created uniform field
        virtual tmp<volScalarField> T() const = 0;

        //- Return the scalar value for isothermal models.
        virtual const dimensionedScalar& Tvalue() const;


        inline const word& modelName() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaEnergyModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

