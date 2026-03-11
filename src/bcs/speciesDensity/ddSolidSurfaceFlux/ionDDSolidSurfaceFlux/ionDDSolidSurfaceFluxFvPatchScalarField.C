/*---------------------------------------------------------------------------*\
  File: ionDDSolidSurfaceFluxFvPatchScalarField.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ionDDSolidSurfaceFluxFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "ionDDSolidSurfaceFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    ionDDSolidSurfaceFluxFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> ionDDSolidSurfaceFluxFvPatchScalarField::calcWallVelocity
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
ionDDSolidSurfaceFluxFvPatchScalarField::
ionDDSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddSolidSurfaceFluxFvPatchScalarField(p, iF),
    enableSurfaceCharging_(true)
{}

// Dictionary Constructor
ionDDSolidSurfaceFluxFvPatchScalarField::
ionDDSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddSolidSurfaceFluxFvPatchScalarField(p, iF, dict),
    enableSurfaceCharging_(dict.lookupOrDefault<bool>
        ("enableSurfaceCharging", true))
{}

// Mapping Constructor
ionDDSolidSurfaceFluxFvPatchScalarField::
ionDDSolidSurfaceFluxFvPatchScalarField
(
    const ionDDSolidSurfaceFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddSolidSurfaceFluxFvPatchScalarField(ptf, p, iF, mapper),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor
ionDDSolidSurfaceFluxFvPatchScalarField::
ionDDSolidSurfaceFluxFvPatchScalarField
(
    const ionDDSolidSurfaceFluxFvPatchScalarField& ptf
)
:
    ddSolidSurfaceFluxFvPatchScalarField(ptf),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor (with new internal field)
ionDDSolidSurfaceFluxFvPatchScalarField::
ionDDSolidSurfaceFluxFvPatchScalarField
(
    const ionDDSolidSurfaceFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddSolidSurfaceFluxFvPatchScalarField(ptf, iF),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ionDDSolidSurfaceFluxFvPatchScalarField::write(Ostream& os) const
{
    ddSolidSurfaceFluxFvPatchScalarField::write(os);   

    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
