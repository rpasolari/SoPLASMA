/*---------------------------------------------------------------------------*\
  File: electromagneticsModel.C
  Part of: SoPLASMA
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

//- Construct from meshes
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
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    mesh_(mesh),
    dielectricMeshes_(dielectricMeshes),
    ePotential_
    (
        IOobject
        (
            "ePotential",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero)
    ),
    Emag_
    (
        IOobject
        (
            "Emag",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mag(E_)
    ),
    phiE_
    (
        IOobject
        (
            "phiE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -fvc::snGrad(ePotential_) * mesh.magSf()
    ),
    reducedE_
    (
        IOobject
        (
            "reducedE",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 4, -3, 0, 0, -1, 0), 0.0)
    ),
    chargeDensity_
    (
        IOobject
        (
            "chargeDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -3, 1, 0, 0, 1, 0), 0.0)
    ),
    surfCharge_
    (
        IOobject
        (
            "surfCharge",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), 0.0)
    ),
    epsilon_
    (
        "epsilon",
        dimensionSet(-1, -3, 4, 0, 0, 2, 0),
        0.0
    ),
    epsilonR_(1.0)
{
    // --- DEBUG PRINTING START ---
    Info<< "\n" << "========================================" << endl;
    Info<< "DEBUG: Checking Object Registries" << endl;
    
    // 1. Check primary mesh registry
    Info<< "Primary Mesh (" << mesh_.name() << ") Registry Contents:" << nl
        << mesh_.thisDb().sortedToc() << nl << endl;

    // 2. Check each dielectric mesh registry
    forAll(dielectricMeshes_, i)
    {
        Info<< "Dielectric Mesh (" << dielectricMeshes_[i].name() << ") Registry Contents:" << nl
            << dielectricMeshes_[i].thisDb().sortedToc() << nl << endl;
    }
    Info<< "========================================\n" << endl;
    // --- DEBUG PRINTING END ---
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<electromagneticsModel> electromagneticsModel::New
(
    const fvMesh& mesh,
    const UPtrList<fvMesh>& dielectricMeshes
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
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
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
        ctorPtr(mesh, dielectricMeshes)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
