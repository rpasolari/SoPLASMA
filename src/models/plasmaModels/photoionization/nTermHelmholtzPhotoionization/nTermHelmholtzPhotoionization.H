/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
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
    Foam::nTermHelmholtzPhotoionization
    
Description

    Implementation of an N-term Helmholtz-equation-based 
    photoionization model for plasma simulations.

    The non-local Zheleznyak integral for the photoionization source is
    approximated by decomposing the absorption kernel into a user-defined
    sum of N exponentials. This yields N independent screened-Poisson 
    (Helmholtz) PDEs solved via the native OpenFOAM Laplacian laplacian():
    
        laplacian(S_j) - sqr(lambda_j * p) * S_j = -A_j * sqr(p) * I

    where:
        S_j      = Photoionization source component for the j-th term
        lambda_j = Fitting parameter lambda for the j-th term [1/(m Pa)]
        A_j      = Fitting parameter A for the j-th term [1/(m Pa2)]
        p        = Total gas pressure [Pa]
        I        = Ionization production rate    

    The total photoionization source term is obtained by summing all N terms:
        S_ph = sum_{j=1}^{N} S_j

    The coefficients A_j and lambda_j are specified by the user obtained via a
    least-squares fit to the generalized photoionization function:
        Phi_0(r) / p = (p * r) * sum_{j=1}^{N} [ A_j * exp(-lambda_j * p * r) ]

    Two unit conventions are supported for the fitting coefficients and the
    gas pressure, selected by the 'coefficientUnits' keyword:

      - SI     : p [Pa],   lambda [1/(m  Pa)],   A [1/(m^2 Pa^2)]
      - TorrCm   : p [Torr], lambda [1/(cm Torr)], A [1/(cm^2 Torr^2)]

    Internally, all quantities are converted to and stored in SI.

Usage
    Example of the model in constant/photoionizationProperties:
    \verbatim
    photoionizationModel    nTermHelmholtzPhotoionization;
    solveInterval           10;         // Default : 1 step

    nTermHelmholtzPhotoionizationCoeffs
    {
        // Unit convention for p, lambda and A: (SI | TorrCm)
        coefficientUnits    TorrCm;

        // Uniform gas pressure (in the chosen unit system)
        p                   760;        // [Torr]

        // Name of the ionization source field in the registry
        I                   I_ion;      // [1/(m3 s)]

        // Number of Helmholtz terms (must match the list size below)
        nTerms              3;

        // Fitting parameters specified as pairs of (lambda A)
        // in the chosen unit system
        fittingParameters
        (
            (0.0553   1.986e-4)
            (0.1460   0.0051)
            (0.8900   0.4886)
        );
    }
    \endverbatim

Reference
    Zheleznyak M B, Mnatsakanyan A Kh, and Sizykh S V (1982)
    Photoionization of nitrogen and oxygen mixtures by radiation
    from a gas discharge. High Temperature Science, 20(3), 357-362.

    Bourdon A, Pasko V P, Liu N Y, Celestin S, Segur P, and Marode E (2007)
    Efficient models for photoionization produced by non-thermal gas
    discharges in air based on radiative transfer and the Helmholtz
    equations. Plasma Sources Science and Technology, 16(3), 656-678.    

SourceFiles
    nTermHelmholtzPhotoionization.C

\*---------------------------------------------------------------------------*/

#ifndef nTermHelmholtzPhotoionization_H
#define nTermHelmholtzPhotoionization_H

#include "PtrList.H"
#include "scalarList.H"
#include "Tuple2.H"
#include "volFields.H"

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class nTermHelmholtzPhotoionization Declaration
\*---------------------------------------------------------------------------*/

class nTermHelmholtzPhotoionization
:
    public photoionizationModel
{
private:

    //- Disallow copy
    nTermHelmholtzPhotoionization
    (
        const nTermHelmholtzPhotoionization&
    ) = delete;

    void operator=(const nTermHelmholtzPhotoionization&) = delete;

    // Private Member Data

        //- Unit convention for the fitting coefficients ("SI" or "TorrCm")
        word coefficientUnits_;

        //- Number of Helmholtz terms (= length of A_ and lambda_)
        label nTerms_;

        //- Constant uniform gas pressure [Pa] (converted to SI internally)
        dimensionedScalar p_;

        //- Name of the ionization production rate field
        word IName_;

        //- Amplitudes A_j in SI units [1/(m2 Pa2)]
        scalarList A_;

        //- Inverse absorption lengths lambda_j in SI units [1/(m Pa)]
        scalarList lambda_;

        //- Partial source fields, one per Helmholtz term
        PtrList<volScalarField> Sph_j_;

public:

    //- Runtime type information
    TypeName("nTermHelmholtzPhotoionization");

    // Constructors

        //- Construct from mesh and dictionary
        nTermHelmholtzPhotoionization(const fvMesh& mesh);

    //- Destructor
    virtual ~nTermHelmholtzPhotoionization() = default;


    // Member Functions

        //- Assemble and solve the N Helmholtz equations, sum into total field
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

