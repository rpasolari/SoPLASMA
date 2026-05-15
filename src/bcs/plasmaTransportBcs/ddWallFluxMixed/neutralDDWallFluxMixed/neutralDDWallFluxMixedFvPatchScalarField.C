/*---------------------------------------------------------------------------*\
  File: neutralDDWallFluxMixedFvPatchScalarField.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::neutralDDWallFluxMixedFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "neutralDDWallFluxMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    neutralDDWallFluxMixedFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> neutralDDWallFluxMixedFvPatchScalarField::calcWallVelocity
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
neutralDDWallFluxMixedFvPatchScalarField::
neutralDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF)
{}

// Dictionary Constructor
neutralDDWallFluxMixedFvPatchScalarField::
neutralDDWallFluxMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddWallFluxMixedFvPatchScalarField(p, iF, dict)
{}

// Mapping Constructor
neutralDDWallFluxMixedFvPatchScalarField::
neutralDDWallFluxMixedFvPatchScalarField
(
    const neutralDDWallFluxMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, p, iF, mapper)
{}

// Copy Constructor
neutralDDWallFluxMixedFvPatchScalarField::
neutralDDWallFluxMixedFvPatchScalarField
(
    const neutralDDWallFluxMixedFvPatchScalarField& ptf
)
:
    ddWallFluxMixedFvPatchScalarField(ptf)
{}

// Copy Constructor (with new internal field)
neutralDDWallFluxMixedFvPatchScalarField::
neutralDDWallFluxMixedFvPatchScalarField
(
    const neutralDDWallFluxMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddWallFluxMixedFvPatchScalarField(ptf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void neutralDDWallFluxMixedFvPatchScalarField::write(Ostream& os) const
{
    ddWallFluxMixedFvPatchScalarField::write(os);   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
