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
    Foam::plasmaPropertyEvaluator

Description:
    Base class for mathematical property evaluation.

    The plasmaPropertyEvaluator class provides a unified interface for 
    mathematical models (interpolation, analytical functions, etc.) that are 
    decoupled from specific physical identities. 

SourceFiles:
    plasmaPropertyEvaluator.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaPropertyEvaluator_H
#define Foam_plasmaPropertyEvaluator_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "plasmaSpecies.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plasmaPropertyEvaluator Declaration
\*---------------------------------------------------------------------------*/

class plasmaPropertyEvaluator
{
protected:

    // Protected Data

        //- The specific math type (e.g., "tabulated1D")
        const word modelType_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the species database
        const plasmaSpecies& species_;

        //- Index of the species this evaluator serves
        const label specieIndex_;

        //- The name of the evaluated field (for prints)
        const word propertyName_;

        //- The dimensionSet of the evaluated field
        const dimensionSet dims_;

public:

    //- Runtime type information
    TypeName("plasmaPropertyEvaluator");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        plasmaPropertyEvaluator,
        dictionary,
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const word& propertyName,
            const dimensionSet& dims
        ),
        (modelType, dict, mesh, species, specieIndex, propertyName, dims)
    );

    // Constructors

        //- Construct from components
        plasmaPropertyEvaluator
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const word& propertyName,
            const dimensionSet& dims
        );

    //- Selector
    static autoPtr<plasmaPropertyEvaluator> New
    (
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh,
        const plasmaSpecies& species,
        const label specieIndex,
        const word& propertyName,
        const dimensionSet& dims
    );

    //- Destructor
    virtual ~plasmaPropertyEvaluator() = default;

    // Member Functions

        // Generate and return a temporary field property
        virtual tmp<volScalarField> evaluate() const = 0;

        //- Return the scalar value for uniform field properties
        virtual dimensionedScalar evaluateScalar() const = 0;

        //- Update an existing field
        virtual void correct(volScalarField& field) const = 0;

        //- Evaluate for a patch field
        virtual void evaluate
        (
            fvPatchScalarField& pField, 
            const label patchi
        ) const = 0;

        //- Returns true if the property is spatially uniform
        virtual bool isUniform() const = 0;

        //- Returns true if the property is constant in time
        virtual bool isConstant() const = 0;

        //- Interface compatibility: Keyword to look for in 'constant' type
        static word coeffKeyword() { return "value"; }

        //- Interface compatibility: Returns the physical model category
        static word propertyType() { return "evaluator"; }

        //- Access to model type string
        word modelType() const { return modelType_; }

        //- Access the name of the field that is evaluated (for prints)
        const word& propertyName() const { return propertyName_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
