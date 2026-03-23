/*---------------------------------------------------------------------------*\
  File: electromagneticsModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::electromagneticsModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "electromagneticsModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(electromagneticsModel, 0);
defineRunTimeSelectionTable(electromagneticsModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from meshes and species
electromagneticsModel::electromagneticsModel
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes,
    const plasmaSpecies* species
)
:
    IOdictionary
    (
        IOobject
        (
            "electromagneticsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    dielectricMeshes_(dielectricMeshes),
    speciesPtr_(species)
{}

//- Construct from meshes only
electromagneticsModel::electromagneticsModel
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
)
:
    IOdictionary
    (
        IOobject
        (
            "electromagneticsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    dielectricMeshes_(dielectricMeshes),
    speciesPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//- Selector for meshes and species
autoPtr<electromagneticsModel> electromagneticsModel::New
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes,
    const plasmaSpecies* species
)
{
    IOdictionary tmpDict
    (
        IOobject
        (
            "electromagneticsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word modelName(tmpDict.get<word>("electromagneticsModel"));

    // Look up the constructor in the table
    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(tmpDict)
            << "Unknown electromagneticsModel type '" << modelName << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<electromagneticsModel>
    (
        ctorPtr(mesh, dielectricMeshes, species)
    );
}

//- Selector for meshes only
autoPtr<electromagneticsModel> electromagneticsModel::New
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
)
{
    return New(mesh, dielectricMeshes, nullptr);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
