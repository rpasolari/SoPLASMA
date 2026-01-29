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
    Foam::fvm::ScharfetterGummel

Description

    Scharfetter-Gummel finite volume matrix operator for drift-diffusion.

    This function builds a coupled fvMatrix for the transport equation:
        div(J) = Source

    Where the total advective-diffusive flux J is defined as:
        J = phiDrift * rho - D * grad(rho)

    This scheme solves the 1D convection-diffusion equation analytically 
    between cell centers (assuming constant E and D locally) to determine 
    the exponential profile of the scalar field. This ensures positivity 
    of the species density even in regions of high electric fields (sheaths) 
    where standard central or upwind schemes may fail or produce oscillations.

    The matrix coefficients are derived from the flux expression:
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
    
\*---------------------------------------------------------------------------*/

#ifndef ScharfetterGummel_H
#define ScharfetterGummel_H

#include "fvMatrix.H"
#include "fvm.H"
#include "fvc.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Function ScharfetterGummel
\*---------------------------------------------------------------------------*/

namespace fvm
{

    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
    (
        const GeometricField<Type, fvPatchField, volMesh>& rho,
        const surfaceScalarField& phiDrift,
        const volScalarField& D
    )
    {
        const fvMesh& mesh = rho.mesh();

        // Initialize the matrix
        tmp<fvMatrix<Type>> tmat
            (
                new fvMatrix<Type>(rho, dimVol/dimTime * rho.dimensions())
            );
        fvMatrix<Type>& mat = tmat.ref();

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
            if (mag(x) < 1e-2) 
            {
                return 1.0 - x/2.0 + sqr(x)/12.0 - pow4(x)/720.0;
            }
            // Limit 2: B(80) = 0 effectively
            if (x > 80)    return 0.0;

            // Limit 3: At x -> -inf, B(x) = -x
            if (x < -80)   return -x;

            // Standard calculation
            return x / (exp(x) - 1.0); 
        };

        // Internal Face loop
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

            mat.diag()[P] += coeffP;
            mat.diag()[N] += coeffN;

            mat.upper()[facei] -= coeffN;
            mat.lower()[facei] -= coeffP;
        }

        // Boundary condition handling
        forAll(rho.boundaryField(), patchi)
        {
            const fvPatchField<Type>& ptf = rho.boundaryField()[patchi];

            scalarField& intCoeffs = mat.internalCoeffs()[patchi];
            scalarField& bdryCoeffs = mat.boundaryCoeffs()[patchi];

            const auto& magSf = mesh.magSf().boundaryField()[patchi];
            const auto& delta = mesh.deltaCoeffs().boundaryField()[patchi];
            const auto& DfPatch = Df.boundaryField()[patchi];
            const auto& phiDriftPatch = phiDrift.boundaryField()[patchi];

            // We only fetch these for physical boundaries.
            // Calling them on coupled patches is unnecessary and can be unsafe.
            const scalarField* alphaPtr = nullptr;
            const scalarField* betaPtr = nullptr;
            tmp<scalarField> tAlpha;
            tmp<scalarField> tBeta;

            if (!ptf.coupled())
            {
                tAlpha = ptf.valueInternalCoeffs(rho);
                tBeta  = ptf.valueBoundaryCoeffs(rho);
                alphaPtr = &tAlpha();
                betaPtr  = &tBeta();
            }

            forAll(ptf, facei)
            {
                // Calculate local properties
                scalar diffCond = DfPatch[facei] * magSf[facei] * delta[facei];
                scalar driftFlux = phiDriftPatch[facei];
                scalar Pe = driftFlux / (diffCond + VSMALL);

                scalar coeffP = diffCond * Bern(-Pe);
                scalar coeffB = diffCond * Bern(Pe);
                if (ptf.coupled())
                {
                    intCoeffs[facei] = coeffP;
                    bdryCoeffs[facei] = coeffB;
                }
                else
                {
                    const scalar alpha = (*alphaPtr)[facei];
                    const scalar beta = (*betaPtr)[facei];

                    intCoeffs[facei] = coeffP - (coeffB * alpha);
                    bdryCoeffs[facei] = coeffB * beta;        
                }
            }
        }

        // Ensure this is okay!!
        mat += fvc::laplacian(Df, rho) - fvm::laplacian(Df, rho);

        return tmat;
    }

    // Overload for diffusivity as a dimensionedScalar
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
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
    tmp<fvMatrix<Type>> ScharfetterGummel
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

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
