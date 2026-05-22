/*---------------------------------------------------------------------------*\
  File: immobileTransportModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::immobileTransportModel.

  Copyright (C) 2026 Rention Pasolari
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
    plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaTransportModel
    (
        modelName, 
        dict, 
        mesh, 
        species, 
        specieIndex
    )
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void immobileTransportModel::correct()
{
    // No tranport coefficients to update
}

tmp<fvScalarMatrix> immobileTransportModel::nEqn() const
{
    return fvm::ddt(species_.numberDensity(specieIndex_));
}

void immobileTransportModel::updateFluxes
(
    const fvScalarMatrix& nEqnMatrix,
    surfaceScalarField&,
    surfaceScalarField&,
    surfaceScalarField&
) const
{
    FatalErrorInFunction
        << "updateFluxes() called for immobile species '"
        << species_.speciesName(specieIndex_) << "'." << nl
        << "Immobile species have no particle flux." << nl
        << "This is a logic error in plasmaTransport." << nl
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
