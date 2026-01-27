/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fvc::ScharfetterGummel

Description

    Scharfetter-Gummel finite volume calculus operator for drift-diffusion.

    This function calculates the explicit advective-diffusive face flux field
    based on the exponential profile assumption:
        J_f = phiDrift * rho - D * grad(rho)

    Unlike the fvm equivalent, this function does not build a matrix.
    It returns the evaluated surfaceScalarField (flux) directly:
        J_f = D/h * [ B(-Pe) * rho_P  -  B(Pe) * rho_N ]

    Where:
      - Pe is the local Peclet number: phiDrift / (D/h)
      - B(x) is the Bernoulli function: x / (exp(x) - 1)
      - rho_P, rho_N are the cell-center values of owner and neighbour.

    For a detailed theoretical derivation, please refer to:
    docs/theory/scharfetter_gummel.md

References
    Scharfetter, D. L., & Gummel, H. K. (1969).
    "Large-Signal Analysis of a Silicon Read Diode Oscillator."
    IEEE Transactions on Electron Devices, 16(1), 64–77.

SourceFiles
    fvcScharfetterGummel.C
    
\*---------------------------------------------------------------------------*/

#ifndef fvcScharfetterGummel_H
#define fvcScharfetterGummel_H

#include "fvMesh.H"
#include "fvc.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Function ScharfetterGummel
\*---------------------------------------------------------------------------*/

namespace fvc
{

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phiDrift,
        const volScalarField& D
    )
    {
        const fvMesh& mesh = rho.mesh();

        // Create the result flux field
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tFlux
        (
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                IOobject
                (
                    "flux(" + rho.name() + ")_SG",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<Type>
                ("flux", phiDrift.dimensions()*rho.dimensions(), 0)
            )
        );

        auto& flux = tFlux.ref();

        const surfaceScalarField& magSf = mesh.magSf();
        const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

        const labelList& owner = mesh.owner();
        const labelList& neighbour = mesh.neighbour();

        // Interpolate diffusivity
        surfaceScalarField Df = fvc::interpolate(D);

        // Bernoulli Function Definition: B(x) = x / exp(x) - 1
        auto Bern = [](scalar x) -> scalar
        {
            // Limit 1: Singularity at x = 0
            if (mag(x) < 1e-4) 
            {
                return 1.0 - x/2.0 + sqr(x)/12.0 - pow4(x)/720.0;
            }
            // Limit 2: B(200) = 0 effectively
            if (x > 200)    return 0.0;

            // Limit 3: At x -> -inf, B(x) = -x
            if (x < -200)   return -x;

            // Standard calculation
            return x / (exp(x) - 1.0); 
        };

        // Internal Face Loop
        for (label facei = 0; facei < owner.size(); ++facei)
        {
            label P = owner[facei];
            label N = neighbour[facei];

            scalar diffCond = Df[facei] * magSf[facei] * deltaCoeffs[facei];
            scalar driftFlux = phiDrift[facei];

            // Local Peclet number
            scalar Pe = driftFlux / (diffCond + VSMALL);

            // Bernoulli weights
            scalar coeffP = diffCond * Bern(-Pe);
            scalar coeffN = diffCond * Bern(Pe);

            // Explicit flux calculation
            flux[facei] = (coeffP * rho[P]) - (coeffN * rho[N]);
        }

        // Boundary Face Loop
        forAll(rho.boundaryField(), patchi)
        {
            const fvPatchField<Type>& rhoBoundary = rho.boundaryField()[patchi];
            tmp<scalarField> tRhoInternal = rhoBoundary.patchInternalField();
            const scalarField& rhoInternal = tRhoInternal();

            auto& fluxPatch = flux.boundaryFieldRef()[patchi];

            const auto& magSf = mesh.magSf().boundaryField()[patchi];
            const auto& delta = mesh.deltaCoeffs().boundaryField()[patchi];
            const auto& DfPatch = Df.boundaryField()[patchi];
            const auto& phiDriftPatch = phiDrift.boundaryField()[patchi];

            forAll(rhoBoundary, facei)
            {
                scalar diffCond = DfPatch[facei] * magSf[facei] *  delta[facei];
                scalar driftFlux = phiDriftPatch[facei];
                scalar Pe = driftFlux / (diffCond + VSMALL);

                scalar coeffP = diffCond * Bern(-Pe);
                scalar coeffB = diffCond * Bern(Pe);

                fluxPatch[facei] = 
                    (coeffP * rhoInternal[facei]) -
                    (coeffB * rhoBoundary[facei]);
            }
        }

        return tFlux;
    }

    // Overload for diffusivity as a dimensionedScalar
    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phiDrift,
        const dimensionedScalar& D
    )
    {
        volScalarField DField
        (
            IOobject
            (
                "D_tmp", 
                rho.mesh().time().timeName(), 
                rho.mesh(), 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ),
            rho.mesh(),
            D
        );    
        
        // Call the core implementation
        return ScharfetterGummel(rho, phiDrift, DField);
    }

    // Overload for Electric Potential and Mobility (mu) instead of flux
    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const volScalarField& ePotential,
        const volScalarField& mu,
        const volScalarField& D,
        const scalar chargeNumber
    )    
    {   
        const fvMesh& mesh = rho.mesh();
        
        surfaceScalarField phiDrift = chargeNumber * fvc::interpolate(mu)
                                    * (-fvc::snGrad(ePotential))
                                    * mesh.magSf();   
  
        // Call the core implementation
        return ScharfetterGummel(rho, phiDrift, D);
    }

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
