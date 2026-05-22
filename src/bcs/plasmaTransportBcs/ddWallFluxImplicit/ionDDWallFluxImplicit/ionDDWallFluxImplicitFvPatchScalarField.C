/*---------------------------------------------------------------------------*\
  File: ionDDWallFluxImplicitFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::ionDDWallFluxImplicitFvPatchScalarField.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

#include "ionDDWallFluxImplicitFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    ionDDWallFluxImplicitFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
ionDDWallFluxImplicitFvPatchScalarField::calcWallThermalVelocity
(
    const dimensionedScalar& m,
    const scalarField& T
) const
{
    // For ions, the thermal velocity is simply u_th/4
    // No drift correction here. Drift is handled separately in
    // valueInternalCoeffs() via the pos()/neg() switch
    return calcThermalVelocity(m, T);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
ionDDWallFluxImplicitFvPatchScalarField::
ionDDWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxImplicitFvPatchScalarField(p, iF),
    enableSurfaceCharging_(false)
{}

// Dictionary Constructor
ionDDWallFluxImplicitFvPatchScalarField::
ionDDWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxImplicitFvPatchScalarField(p, iF, dict),
    enableSurfaceCharging_
            (dict.lookupOrDefault<bool>("enableSurfaceCharging", false))
{}

// Mapping Constructor
ionDDWallFluxImplicitFvPatchScalarField::
ionDDWallFluxImplicitFvPatchScalarField
(
    const ionDDWallFluxImplicitFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf, p, iF, mapper),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor
ionDDWallFluxImplicitFvPatchScalarField::
ionDDWallFluxImplicitFvPatchScalarField
(
    const ionDDWallFluxImplicitFvPatchScalarField& ptf
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// Copy Constructor (with new internal field)
ionDDWallFluxImplicitFvPatchScalarField::
ionDDWallFluxImplicitFvPatchScalarField
(
    const ionDDWallFluxImplicitFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf, iF),
    enableSurfaceCharging_(ptf.enableSurfaceCharging_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ionDDWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxImplicitFvPatchScalarField::write(os);   

    os.writeEntry("enableSurfaceCharging", enableSurfaceCharging_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
