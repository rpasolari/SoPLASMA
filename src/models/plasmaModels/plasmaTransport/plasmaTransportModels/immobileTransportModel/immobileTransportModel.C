/*---------------------------------------------------------------------------*\
  File: immobileTransportModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::immobileTransportModel.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "immobileTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(immobileTransportModel, 0);
addToRunTimeSelectionTable
(
    plasmaTransportModel,
    immobileTransportModel, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immobileTransportModel::immobileTransportModel
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
    plasmaTransportModel
    (
        modelName, 
        dict, 
        mesh, 
        species, 
        specieIndex, 
        E,
        phiE
    )
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

tmp<fvScalarMatrix> immobileTransportModel::nEqn() const
{
    return tmp<fvScalarMatrix>(fvm::ddt(species_.numberDensity(specieIndex_)));
}

void immobileTransportModel::updateParticleFlux(surfaceScalarField& flux) const
{
    FatalErrorInFunction
        << "Attempted to update particle flux for immobile species '"
        << species_.speciesNames()[specieIndex_] << "'." << nl
        << "This indicates a logic error in the plasmaTransport class."
        << abort(FatalError);
}

tmp<volScalarField> immobileTransportModel::electricalConductivity() const
{
    return tmp<volScalarField>::New
    (
        IOobject("sigma0", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("0", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0)
    );
}

tmp<volScalarField> immobileTransportModel::diffusiveChargeSource() const
{
    return tmp<volScalarField>::New
    (
        IOobject("diffSrc0", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0.0)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
