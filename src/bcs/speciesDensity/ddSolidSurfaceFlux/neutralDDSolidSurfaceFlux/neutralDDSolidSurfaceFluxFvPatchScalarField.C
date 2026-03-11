/*---------------------------------------------------------------------------*\
  File: neutralDDSolidSurfaceFluxFvPatchScalarField.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::neutralDDSolidSurfaceFluxFvPatchScalarField.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "neutralDDSolidSurfaceFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Register the class to the runtime selection table
makePatchTypeField
(
    fvPatchScalarField,
    neutralDDSolidSurfaceFluxFvPatchScalarField
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> neutralDDSolidSurfaceFluxFvPatchScalarField::calcWallVelocity
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
neutralDDSolidSurfaceFluxFvPatchScalarField::
neutralDDSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddSolidSurfaceFluxFvPatchScalarField(p, iF)
{}

// Dictionary Constructor
neutralDDSolidSurfaceFluxFvPatchScalarField::
neutralDDSolidSurfaceFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    ddSolidSurfaceFluxFvPatchScalarField(p, iF, dict)
{}

// Mapping Constructor
neutralDDSolidSurfaceFluxFvPatchScalarField::
neutralDDSolidSurfaceFluxFvPatchScalarField
(
    const neutralDDSolidSurfaceFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ddSolidSurfaceFluxFvPatchScalarField(ptf, p, iF, mapper)
{}

// Copy Constructor
neutralDDSolidSurfaceFluxFvPatchScalarField::
neutralDDSolidSurfaceFluxFvPatchScalarField
(
    const neutralDDSolidSurfaceFluxFvPatchScalarField& ptf
)
:
    ddSolidSurfaceFluxFvPatchScalarField(ptf)
{}

// Copy Constructor (with new internal field)
neutralDDSolidSurfaceFluxFvPatchScalarField::
neutralDDSolidSurfaceFluxFvPatchScalarField
(
    const neutralDDSolidSurfaceFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    ddSolidSurfaceFluxFvPatchScalarField(ptf, iF)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void neutralDDSolidSurfaceFluxFvPatchScalarField::write(Ostream& os) const
{
    ddSolidSurfaceFluxFvPatchScalarField::write(os);   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
