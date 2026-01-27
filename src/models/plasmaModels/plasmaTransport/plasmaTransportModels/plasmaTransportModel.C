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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransportModel::plasmaTransportModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E,
    const volScalarField& ePotential
)
:
    modelName_(modelName),
    mesh_(mesh),
    species_(species),
    dict_(dict),
    specieIndex_(specieIndex),
    E_(E),
    ePotential_(ePotential),
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
    const volVectorField& E,
    const volScalarField& ePotential
)
{
    if (!zeroVectorFieldPtr_.valid())
    {
        zeroVectorFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "zeroDriftVelocity",
                    mesh.time().constant(), 
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector
                (
                    "zero", 
                    dimensionSet(0, 1, -1, 0, 0, 0, 0),
                    vector::zero
                )
            )
        );
    }

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
        ctorPtr(modelName, dict, mesh, species, specieIndex, E, ePotential)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
