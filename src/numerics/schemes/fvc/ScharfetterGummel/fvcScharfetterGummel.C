/*---------------------------------------------------------------------------*\
  File: fvcScharfetterGummel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Explicit template instantiation for the ScharfetterGummel finite 
    volume calculus.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
    See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvcScharfetterGummel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fvc
{
    // Explicit instantiation for scalar fields

        //- Core implementation
        template tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>
        ScharfetterGummel
        (
            const GeometricField<scalar, fvPatchField, volMesh>&,
            const surfaceScalarField&,
            const volScalarField&
        );

        // Overload for diffusivity as a dimensionedScalar
        template tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>
        ScharfetterGummel
        (
            const GeometricField<scalar, fvPatchField, volMesh>&,
            const surfaceScalarField&,
            const dimensionedScalar&
        );        

        // Overload for Electric Potential and Mobility (mu) instead of flux
        template tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>>
        ScharfetterGummel
        (
            const GeometricField<scalar, fvPatchField, volMesh>&,
            const volScalarField&,
            const volScalarField&,
            const volScalarField&,
            const scalar
        );

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
