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
        J_f = phi * rho - D * grad(rho)

    Unlike the fvm equivalent, this function does not build a matrix.
    It returns the evaluated surfaceScalarField (flux) directly:
        J_f = D/h * [ B(-Pe) * rho_P  -  B(Pe) * rho_N ]

    Where:
      - Pe is the local Peclet number: phi / (D/h)
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
        const surfaceScalarField& phi,
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
                ("flux", phi.dimensions()*rho.dimensions(), Zero)
            )
        );

        auto& flux = tFlux.ref();

        const surfaceScalarField& magSf = mesh.magSf();
        const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

        const labelList& owner = mesh.owner();
        const labelList& neighbour = mesh.neighbour();

        surfaceScalarField Df = fvc::interpolate(D);

        // Bernoulli Function Definition: B(x) = x / exp(x) - 1
        auto Bern = [](scalar x) -> scalar
        {
            const scalar ax = mag(x);
            if (ax < 1e-4)
            {
                return 1.0 - 0.5*x + (x*x)/12.0 - pow4(x)/720.0;;
            }
            if (x > 200.0)  return 0.0;
            if (x < -200.0) return -x;
            return x / (Foam::exp(x) - 1.0);
        };

        const scalarField& convFlux = phi.primitiveField();
        const scalarField diffCond = 
            Df.primitiveField() * magSf.primitiveField()
          * deltaCoeffs.primitiveField();

        // Internal Face Loop
        for (label facei = 0; facei < owner.size(); ++facei)
        {
            label P = owner[facei];
            label N = neighbour[facei];

            scalar Pe = convFlux[facei] / (diffCond[facei] + VSMALL);

            scalar coeffP = diffCond[facei] * Bern(-Pe);
            scalar coeffN = diffCond[facei] * Bern(Pe);

            flux[facei] = (coeffP * rho[P]) - (coeffN * rho[N]);
        }

        // Boundary Face Loop
        forAll(rho.boundaryField(), patchi)
        {
            const fvPatchField<Type>& ptf = rho.boundaryField()[patchi];
            
            // Get internal field values for the patch
            tmp<Field<Type>> tRhoInternal = ptf.patchInternalField();
            const Field<Type>& rhoInternal = tRhoInternal();

            auto& fluxPatch = flux.boundaryFieldRef()[patchi];

            const scalarField pDiffCond = Df.boundaryField()[patchi] 
                * mesh.magSf().boundaryField()[patchi] 
                * mesh.deltaCoeffs().boundaryField()[patchi];
                
            const scalarField& pPhi = phi.boundaryField()[patchi];

            forAll(ptf, facei)
            {
                scalar Pe = pPhi[facei] / (pDiffCond[facei] + VSMALL);

                scalar coeffP = pDiffCond[facei] * Bern(-Pe);
                scalar coeffB = pDiffCond[facei] * Bern(Pe);

                // For FVC, we use the actual boundary values (ptf[facei])
                fluxPatch[facei] = 
                    (coeffP * rhoInternal[facei]) -
                    (coeffB * ptf[facei]);
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
        const surfaceScalarField& phi,
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
        return ScharfetterGummel(rho, phi, DField);
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
        
        surfaceScalarField phi = chargeNumber * fvc::interpolate(mu)
                                    * (-fvc::snGrad(ePotential))
                                    * mesh.magSf();   
  
        // Call the core implementation
        return ScharfetterGummel(rho, phi, D);
    }

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
