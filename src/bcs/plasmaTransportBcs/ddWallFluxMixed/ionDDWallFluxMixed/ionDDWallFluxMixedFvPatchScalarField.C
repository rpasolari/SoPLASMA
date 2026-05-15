/*---------------------------------------------------------------------------*\
  File: ionDDWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ionDDWallFluxMixedFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "ionDDWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    ionDDWallFluxMixedFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> ionDDWallFluxMixedFvPatchScalarField::calcWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal,
    const scalar Z
) const
{
    // This is mathematically identical to uTh - Z*uD + max(0, Z*uD)
    // but uses half the operations.
    return calcThermalVelocity(m, T) + max(0.0, -Z * uDriftNormal);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
ionDDWallFluxMixedFvPatchScalarField::
ionDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF),
    enableSurfaceCharging_(true)
{}

// Dictionary Constructor
ionDDWallFluxMixedFvPatchScalarField::
ionDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF, dict),
    enableSurfaceCharging_(dict.lookupOrDefault<bool>
        ("enableSurfaceCharging", true))
{}

// Mapping Constructor
ionDDWallFluxMixedFvPatchScalarField::
ionDDWallFluxMixedFvPatchScalarField
(
    const ionDDWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, p, iF, mapper),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor
ionDDWallFluxMixedFvPatchScalarField::
ionDDWallFluxMixedFvPatchScalarField
(
    const ionDDWallFluxMixedFvPatchScalarField& ptf
)
:
    ddWallFluxMixedFvPatchScalarField(ptf),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor (with new internal field)
ionDDWallFluxMixedFvPatchScalarField::
ionDDWallFluxMixedFvPatchScalarField
(
    const ionDDWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, iF),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ionDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxMixedFvPatchScalarField::write(os);   

    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
