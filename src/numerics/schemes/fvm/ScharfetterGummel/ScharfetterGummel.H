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
        J = phi * rho - D * grad(rho)

    This scheme solves the 1D convection-diffusion equation analytically 
    between cell centers (assuming constant E and D locally) to determine 
    the exponential profile of the scalar field. This ensures positivity 
    of the species density even in regions of high electric fields (sheaths) 
    where standard central or upwind schemes may fail or produce oscillations.

    The matrix coefficients are derived from the flux expression:
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
        const surfaceScalarField& phi, 
        const volScalarField& D
    )
    {
        const fvMesh& mesh = rho.mesh();

        tmp<fvMatrix<Type>> tfvm
        (
            new fvMatrix<Type>
            (
                rho,
                phi.dimensions() * rho.dimensions()
            )
        );
        fvMatrix<Type>& fvm = tfvm.ref();
    
        surfaceScalarField Df = fvc::interpolate(D);

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
            Df.primitiveField() * mesh.magSf().primitiveField()
          * mesh.deltaCoeffs().primitiveField();

        const scalarField Pe = convFlux / (diffCond + VSMALL);

        scalarField coeffP(diffCond.size());
        scalarField coeffN(diffCond.size());

        forAll(Pe, facei)
        {
            coeffP[facei] = diffCond[facei] * Bern(-Pe[facei]);
            coeffN[facei] = diffCond[facei] * Bern(Pe[facei]);
        }

        fvm.lower() = -coeffP; 
        fvm.upper() = -coeffN;
        fvm.negSumDiag();

        // Boundaries
        forAll(rho.boundaryField(), patchi)
        {
            const fvPatchField<scalar>& ptf = rho.boundaryField()[patchi];

            scalarField diffCond = Df.boundaryField()[patchi] 
                * mesh.magSf().boundaryField()[patchi] 
                * mesh.deltaCoeffs().boundaryField()[patchi];
                
            scalarField Pe = phi.boundaryField()[patchi] 
                / (diffCond + VSMALL);

            scalarField coeffP(Pe.size());
            scalarField coeffB(Pe.size());  

            forAll(ptf, facei)
            {
                coeffP[facei] = diffCond[facei] * Bern(-Pe[facei]);
                coeffB[facei] = diffCond[facei] * Bern(Pe[facei]);
            }

            if (ptf.coupled())
            {
                fvm.internalCoeffs()[patchi] = coeffP;
                fvm.boundaryCoeffs()[patchi] = coeffB;
            }
            else
            {
                tmp<scalarField> tAlpha = ptf.valueInternalCoeffs(rho);
                tmp<scalarField> tBeta = ptf.valueBoundaryCoeffs(rho);

                fvm.internalCoeffs()[patchi] = coeffP - (coeffB * tAlpha());
                fvm.boundaryCoeffs()[patchi] = coeffB * tBeta();
            }
        }

        return tfvm;
    }

    // Overload for diffusivity as a dimensionedScalar
    template<class Type>
    tmp<fvMatrix<Type>> ScharfetterGummel
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
        
        surfaceScalarField phi = chargeNumber * fvc::interpolate(mu)
                                    * (-fvc::snGrad(ePotential))
                                    * mesh.magSf();   
  
        // Call the core implementation
        return ScharfetterGummel(rho, phi, D);
    }

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
