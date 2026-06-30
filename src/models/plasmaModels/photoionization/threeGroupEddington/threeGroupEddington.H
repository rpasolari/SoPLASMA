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
    Foam::threeGroupEddington
    
Description

    Implementation of the three-group Eddington (first-order spherical
    harmonics, P1) approximation to the radiative transfer equation for
    photoionization source calculation in plasma simulations.

    Unlike the Helmholtz model, which directly fits the Zheleznyak g(R)/p
    kernel, the Eddington model is derived on physical grounds from the
    radiative transfer equation itself. For each spectral group j, the
    isotropic part of the photon distribution function Psi*_j satisfies a
    screened-Poisson (Helmholtz-type) PDE sourced directly by the
    ionization production rate I, solved via the native OpenFOAM
    laplacian():

        laplacian(Psi*_j) - 3 * sqr(lambda_j * p) * Psi*_j
            = -(3 * lambda_j * p / c) * I

    where:
        Psi*_j   = xi-weighted isotropic photon distribution component
                   for the j-th spectral group [1/m2]
        lambda_j = Fitting parameter lambda for the j-th group [1/(m Pa)]
        A_j      = Fitting parameter A for the j-th group [1/(m Pa)]
        p        = Total gas pressure [Pa]
        c        = Speed of light [m/s]
        I        = Ionization production rate [1/(m3 s)]

    The photoionization source contribution of each group is recovered from
    the solved field as:

        S_j  = A_j * p * c * Psi*_j

    and the total photoionization source term is obtained by summing the
    three groups:

        S_ph = sum_{j=1}^{3} S_j

    The coefficients A_j and lambda_j correspond to the three-group fit of
    the generalized photoionization function:
        Phi_0(r) / p = sum_{j=1}^{3} [ A_j * exp(-lambda_j * p * r) ]

    Note the difference in form (and resulting coefficient units) compared
    to the Helmholtz model: here the fit is applied directly to g(R)/p
    without the additional (p * r) prefactor, since the prefactor arises
    naturally from the P1 closure of the radiative transfer equation rather
    than from an ad hoc kernel approximation.

    Boundary conditions
        Supported values of 'boundaryConditionType':

          - marshak : Marshak (Robin) condition derived from radiative
                      transfer theory, applied independently to each
                      Psi*_j:

                          grad(Psi*_j) & n_S = -(3/2) * lambda_j * p * Psi*_j

          - zero    : Zero Dirichlet boundary conditions for all Psi*_j.

        Any other value triggers a FatalError.

    Note
        This class solves the PDE for whatever A_j, lambda_j, p, and I are
        provided; it makes no assumption about the gas. For the
        Zheleznyak/Bourdon air parametrization, p is the O2 partial
        pressure (not the total), and I must already include the quenching
        factor pq/(p+pq) and the xi*nu_u/nu_i ratio. The user is
        responsible for ensuring A_j, lambda_j, p, and I are consistent
        for the gas of interest.

    Two unit conventions are supported for the fitting coefficients and the
    gas pressure, selected by the 'coefficientUnits' keyword:

      - SI     : p [Pa],   lambda [1/(m  Pa)],   A [1/(m Pa)]
      - TorrCm : p [Torr], lambda [1/(cm Torr)], A [1/(cm Torr)]

    Internally, all quantities are converted to and stored in SI.

Usage
    Example of the model in constant/photoionizationProperties:
    \verbatim
    photoionizationModel    threeGroupEddington;
    solveInterval           10;         // Default : 1 step

    threeGroupEddingtonCoeffs
    {
        // Unit convention for p, lambda and A: (SI | TorrCm)
        coefficientUnits        TorrCm;

        // Boundary condition type: (zero | marshak)
        boundaryConditionType   marshak;

        // Uniform gas pressure (in the chosen unit system)
        p                       150;        // [Torr]

        // Name of the ionization source field in the registry
        I                       I_ion;      // [1/(m3 s)]

        // Fitting parameters specified as pairs of (lambda A)
        // in the chosen unit system, for the three spectral groups
        fittingParameters
        (
            (0.0447   0.0067)
            (0.1121   0.0346)
            (0.5994   0.3059)
        );
    }
    \endverbatim

Reference
    Zheleznyak M B, Mnatsakanyan A Kh, and Sizykh S V (1982)
    Photoionization of nitrogen and oxygen mixtures by radiation
    from a gas discharge. High Temperature Science, 20(3), 357-362.

    Segur P, Bourdon A, Marode E, Bessieres D, and Paillol J H (2006)
    The use of an improved Eddington approximation to facilitate the
    calculation of photoionization in streamer discharges.
    Plasma Sources Science and Technology, 15(4), 648-660.

    Bourdon A, Pasko V P, Liu N Y, Celestin S, Segur P, and Marode E (2007)
    Efficient models for photoionization produced by non-thermal gas
    discharges in air based on radiative transfer and the Helmholtz
    equations. Plasma Sources Science and Technology, 16(3), 656-678.

SourceFiles
    threeGroupEddington.C

\*---------------------------------------------------------------------------*/

#ifndef threeGroupEddington_H
#define threeGroupEddington_H

#include "PtrList.H"
#include "scalarList.H"
#include "Tuple2.H"
#include "volFields.H"

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class threeGroupEddington Declaration
\*---------------------------------------------------------------------------*/

class threeGroupEddington
:
    public photoionizationModel
{
private:

    //- Disallow copy
    threeGroupEddington
    (
        const threeGroupEddington&
    ) = delete;

    void operator=(const threeGroupEddington&) = delete;

    // Private Member Data

        //- Unit convention for the fitting coefficients ("SI" or "TorrCm")
        word coefficientUnits_;

        //- Boundary condition type
        word boundaryConditionType_;

        //- Number of Eddington groups (fixed = 3)
        static const label nGroups_ = 3;

        //- Constant uniform gas pressure [Pa] (converted to SI internally)
        dimensionedScalar p_;

        //- Name of the ionization production rate field
        word IName_;

        //- Amplitudes A_j in SI units [1/(m Pa)]
        scalarList A_;

        //- Inverse absorption lengths lambda_j in SI units [1/(m Pa)]
        scalarList lambda_;

        //- Speed of light [m/s]
        dimensionedScalar c_;

        //- Solved photon distribution components Psi*_j, one per group
        PtrList<volScalarField> PsiStar_j_;

        //- Partial source fields, one per Helmholtz term
        PtrList<volScalarField> Sph_j_;

public:

    //- Runtime type information
    TypeName("threeGroupEddington");

    // Constructors

        //- Construct from mesh and dictionary
        threeGroupEddington(const fvMesh& mesh);

    //- Destructor
    virtual ~threeGroupEddington() = default;


    // Member Functions

        //- Assemble and solve the 3 Eddington equations, sum into total field
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

