/*---------------------------------------------------------------------------*\
  File: neutralDDWallFluxImplicitFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::neutralDDWallFluxImplicitFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "neutralDDWallFluxImplicitFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    neutralDDWallFluxImplicitFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> neutralDDWallFluxImplicitFvPatchScalarField::calcWallVelocity
(
    const dimensionedScalar& m,
    const scalarField& T,
    const scalarField& uDriftNormal,
    const scalar Z
) const
{
    return calcThermalVelocity(m, T);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Standard Constructor
neutralDDWallFluxImplicitFvPatchScalarField::
neutralDDWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxImplicitFvPatchScalarField(p, iF)
{}

// Dictionary Constructor
neutralDDWallFluxImplicitFvPatchScalarField::
neutralDDWallFluxImplicitFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxImplicitFvPatchScalarField(p, iF, dict)
{}

// Mapping Constructor
neutralDDWallFluxImplicitFvPatchScalarField::
neutralDDWallFluxImplicitFvPatchScalarField
(
    const neutralDDWallFluxImplicitFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf, p, iF, mapper)
{}

// Copy Constructor
neutralDDWallFluxImplicitFvPatchScalarField::
neutralDDWallFluxImplicitFvPatchScalarField
(
    const neutralDDWallFluxImplicitFvPatchScalarField& ptf
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf)
{}

// Copy Constructor (with new internal field)
neutralDDWallFluxImplicitFvPatchScalarField::
neutralDDWallFluxImplicitFvPatchScalarField
(
    const neutralDDWallFluxImplicitFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxImplicitFvPatchScalarField(ptf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void neutralDDWallFluxImplicitFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxImplicitFvPatchScalarField::write(os);   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
