/*---------------------------------------------------------------------------*\
  File: plasmaTransportModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTransportModel.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaTransportModel, 0);
defineRunTimeSelectionTable(plasmaTransportModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransportModel::plasmaTransportModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E,
    const surfaceScalarField& phiE
)
:
    modelName_(modelName),
    mesh_(mesh),
    species_(species),
    dict_(dict),
    specieIndex_(specieIndex),
    E_(E),
    phiE_(phiE)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaTransportModel> plasmaTransportModel::New
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E,
    const surfaceScalarField& phiE
)
{
    // Lookup constructor using function-call operator
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        const word& sName = species.speciesNames()[specieIndex];

        FatalIOErrorInFunction(dict)
            << "Species '" << sName << "': "
            << "Unknown plasmaTransportModel type '" << modelName << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    // Construct and return the model
    return autoPtr<plasmaTransportModel>
    (
        ctorPtr(modelName, dict, mesh, species, specieIndex, E, phiE)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
