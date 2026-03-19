/*---------------------------------------------------------------------------*\
  File: plasmaEnergyModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaEnergyModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaEnergyModel, 0);
defineRunTimeSelectionTable(plasmaEnergyModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaEnergyModel::plasmaEnergyModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
:
    modelName_(modelName),
    mesh_(mesh),
    species_(species),
    dict_(dict),
    specieIndex_(specieIndex),
    E_(E)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaEnergyModel> plasmaEnergyModel::New
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
{
    // Lookup constructor using function-call operator
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        const word& sName = species.speciesNames()[specieIndex];
        
        FatalIOErrorInFunction(dict)
            << "Species '" << sName << "': "
            << "Unknown plasmaEnergyModel type '" << modelName << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    // Construct and return the model
    return autoPtr<plasmaEnergyModel>
    (
        ctorPtr(modelName, dict, mesh, species, specieIndex, E)
    );
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

const dimensionedScalar& plasmaEnergyModel::Tvalue() const
{
    const word& sName = species_.speciesNames()[specieIndex_];

    // Default behavior: If a derived class is field-based and doesn't 
    // override this, this function throws an error.
    FatalErrorInFunction
        << "Requested Tvalue() for species '" << sName 
        << "' from a non-isothermal energy model type (" << modelName_ << ")."
        << nl << "This species uses a spatial temperature field. "
        << "Please update the solver logic to use the .T() accessor instead."
        << abort(FatalError);

    // Never reached
    static const dimensionedScalar T_dummy("T_dummy", dimTemperature, 0.0);
    return T_dummy;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
