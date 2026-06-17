/*---------------------------------------------------------------------------*\
  File: immobile.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::immobile.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "immobile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(immobile, 0);
addToRunTimeSelectionTable
(
    plasmaTransportModel,
    immobile, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immobile::immobile
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

void immobile::correct()
{
    // No tranport coefficients to update
}

tmp<fvScalarMatrix> immobile::nEqn() const
{
    return fvm::ddt(species_.numberDensity(specieIndex_));
}

void immobile::updateFluxes
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

tmp<volScalarField> immobile::electricalConductivity() const
{
    FatalErrorInFunction
        << "electricalConductivity() called for immobile species '"
        << species_.speciesName(specieIndex_) << "'." << nl
        << "Immobile species have zero conductivity and should be "
        << "skipped by the caller (use mobileSpeciesIDs/chargedSpeciesIDs)." << nl
        << abort(FatalError);
    return nullptr; 
}

tmp<volScalarField> immobile::diffusiveChargeSource() const
{
    FatalErrorInFunction
        << "diffusiveChargeSource() called for immobile species '"
        << species_.speciesName(specieIndex_) << "'." << nl
        << "Immobile species have zero diffusive charge source and should be "
        << "skipped by the caller (use mobileSpeciesIDs/chargedSpeciesIDs)." << nl
        << abort(FatalError);
    return nullptr; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
