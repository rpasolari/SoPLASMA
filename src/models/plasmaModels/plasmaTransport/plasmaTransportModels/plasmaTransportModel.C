/*---------------------------------------------------------------------------*\
  File: plasmaTransportModel.C
  Part of: foamPlasmaToolkit
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

autoPtr<volVectorField> plasmaTransportModel::zeroVectorFieldPtr_(nullptr);
autoPtr<surfaceScalarField> plasmaTransportModel::zeroSurfaceFieldPtr_(nullptr);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaTransportModel::checkSharedFields(const fvMesh& mesh)
{
    if (!zeroSurfaceFieldPtr_.valid())
    {
        zeroSurfaceFieldPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject("phi0", mesh.time().constant(), mesh),
                mesh,
                dimensionedScalar
                (
                    "0", 
                    dimensionSet(0, 3, -1, 0, 0, 0, 0), 
                    0.0
                )
            )
        );
        
        zeroVectorFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject("v0", mesh.time().constant(), mesh),
                mesh,
                dimensionedVector
                (
                    "0", 
                    dimensionSet(0, 1, -1, 0, 0, 0, 0), 
                    vector::zero
                )
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransportModel::plasmaTransportModel
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
    E_(E),
    fluxScheme_("standard")
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaTransportModel> plasmaTransportModel::New
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
{
    checkSharedFields(mesh);

    // Lookup constructor using function-call operator
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown plasmaTransportModel type '" << modelName << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    // Construct and return the model
    return autoPtr<plasmaTransportModel>
    (
        ctorPtr(modelName, dict, mesh, species, specieIndex, E)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
