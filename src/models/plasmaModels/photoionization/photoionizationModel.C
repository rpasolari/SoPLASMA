/*---------------------------------------------------------------------------*\
  File: photoionizationModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::photoionizationModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(photoionizationModel, 0);
defineRunTimeSelectionTable(photoionizationModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from mesh
photoionizationModel::photoionizationModel(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "photoionizationProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    Sph_
    (
        IOobject
        (
            "Sph",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    solveInterval_(getOrDefault<label>("solveInterval", 1))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<photoionizationModel> photoionizationModel::New
(
  const fvMesh& mesh
)
{
    IOdictionary tmpDict
    (
        IOobject
        (
            "photoionizationProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word modelName(tmpDict.get<word>("photoionizationModel"));

    // Look up the constructor in the table
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(tmpDict)
            << "Unknown photoionizationModel type '" << modelName << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<photoionizationModel>(ctorPtr(mesh));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
