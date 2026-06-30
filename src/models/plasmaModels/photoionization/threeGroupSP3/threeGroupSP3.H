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
    Foam::threeGroupSP3

Description

    Implementation of the three-group SP3 (third-order spherical harmonics,
    improved Eddington) approximation to the radiative transfer equation for
    photoionization source calculation in plasma simulations.

    The SP3 approximation improves on the first-order Eddington model when
    the absorption coefficient of the gas is small and the gradient of the
    source term is large. For each spectral group j, the isotropic part of
    the photon distribution function is obtained by combining two auxiliary
    fields phi1_j and phi2_j, each satisfying its own screened-Poisson
    (Helmholtz-type) PDE sourced directly by the ionization production rate
    I, solved via the native OpenFOAM laplacian():

        laplacian(phi1_j) - (lambda_j * p)^2 / kappa1^2 * phi1_j
            = -(lambda_j * p) / (kappa1^2 * c) * I

        laplacian(phi2_j) - (lambda_j * p)^2 / kappa2^2 * phi2_j
            = -(lambda_j * p) / (kappa2^2 * c) * I

    where:
        phi1_j, phi2_j = auxiliary xi-weighted photon distribution
                         components for the j-th spectral group [1/m3]
        lambda_j       = Fitting parameter lambda for the j-th group
                         [1/(m Pa)]
        A_j            = Fitting parameter A for the j-th group [1/(m Pa)]
        p              = Total gas pressure [Pa]
        c              = Speed of light [m/s]
        I              = Ionization production rate [1/(m3 s)]
        kappa1, kappa2 = SP3 closure constants (dimensionless):
                             kappa1^2 = 3/7 - (2/7)*sqrt(6/5)
                             kappa2^2 = 3/7 + (2/7)*sqrt(6/5)

    The two auxiliary fields are combined via the SP3 closure to recover the
    photoionization source contribution of each group:

        Psi_j = (gamma2 * phi1_j - gamma1 * phi2_j) / (gamma2 - gamma1)

        S_j   = A_j * p * c * Psi_j

    where gamma1, gamma2 are the SP3 closure weights (dimensionless):

        gamma_n = (5/7) * [1 + (-1)^n * (3/5) * sqrt(6/5)]

    and the total photoionization source term is obtained by summing the
    three groups:

        S_ph = sum_{j=1}^{3} S_j

    The coefficients A_j and lambda_j correspond to the same three-group fit
    used by the Eddington model (Bourdon et al. 2007, table 3) of the
    generalized photoionization function:
        Phi_0(r) / p = sum_{j=1}^{3} [ A_j * exp(-lambda_j * p * r) ]

    As demonstrated in Bourdon et al. (2007), the SP3 model can be shown to
    be mathematically equivalent to a six-exponential Helmholtz model, and
    generally provides better agreement with the classical Zheleznyak
    integral model than either the three-exponential Helmholtz or the
    three-group Eddington model, particularly close to domain boundaries
    and at large pO2R values.

    Boundary conditions
        Supported values of 'boundaryConditionType':

          - larsen : Coupled Robin BC from radiative transfer theory
                     (Liu et al. 2007), applied jointly to phi1_j and
                     phi2_j:

                         grad(phi1_j) & n_S = -lambda_j*p*( alpha1*phi1_j
                                                         + beta2 *phi2_j )

                         grad(phi2_j) & n_S = -lambda_j*p*( alpha2*phi2_j
                                                         + beta1 *phi1_j )

                     alpha >> beta, so each equation is solved
                     independently: alpha via mixed BC valueFraction,
                     beta as a lagged refGrad, Picard-iterated per
                     correct() call.

          - zero   : Zero Dirichlet boundary conditions for both phi1_j
                     and phi2_j.

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
    photoionizationModel    3GroupSP3;
    solveInterval           10;         // Default : 1 step

    3GroupSP3Coeffs
    {
        // Unit convention for p, lambda and A: (SI | TorrCm)
        coefficientUnits        TorrCm;

        // Boundary condition type: (larsen | zero)
        boundaryConditionType   larsen;

        // Larsen BC inner Picard iteration controls (optional, defaults
        // shown). Ignored for 'zero' BC.
        larsenMaxIters          10;
        larsenTol               1e-6;

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

    Liu N Y, Celestin S, Bourdon A, Pasko V P, Segur P, and Marode E (2007)
    Application of photoionization models based on radiative transfer and
    the Helmholtz equations to studies of streamers in weak electric
    fields. Appl. Phys. Lett., 91 (21), 211501.

    Larsen E W, Thommes G, Klar A, Seaid M, and Gotz T (2002)
    Simplified PN approximations to the equations of radiative heat
    transfer and applications. J. Comput. Phys., 183 (2), 652-675.

SourceFiles
    threeGroupSP3.C

\*---------------------------------------------------------------------------*/

#ifndef threeGroupSP3_H
#define threeGroupSP3_H

#include "PtrList.H"
#include "scalarList.H"
#include "Tuple2.H"
#include "volFields.H"

#include "photoionizationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class threeGroupSP3 Declaration
\*---------------------------------------------------------------------------*/

class threeGroupSP3
:
    public photoionizationModel
{
private:

    //- Disallow copy
    threeGroupSP3(const threeGroupSP3&) = delete;

    void operator=(const threeGroupSP3&) = delete;

    // Private Member Data

        //- Unit convention for the fitting coefficients ("SI" or "TorrCm")
        word coefficientUnits_;

        //- Boundary condition type
        word boundaryConditionType_;

        //- Number of SP3 groups (fixed = 3)
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

        //- SP3 closure constants kappa1^2, kappa2^2 (dimensionless)
        //  kappa1^2 = 3/7 - (2/7)*sqrt(6/5)
        //  kappa2^2 = 3/7 + (2/7)*sqrt(6/5)
        static const scalar kappa1Sqr_;
        static const scalar kappa2Sqr_;

        //- SP3 closure weights gamma1, gamma2 (dimensionless)
        //  gamma_n = (5/7)*[1 + (-1)^n * (3/5)*sqrt(6/5)]
        static const scalar gamma1_;
        static const scalar gamma2_;

        //- Larsen BC coefficients alpha1, alpha2 (dimensionless)
        //  alpha_{1,2} = (5/96) * (34 +- 11*sqrt(6/5))
        static const scalar alpha1_;
        static const scalar alpha2_;

        //- Larsen BC coefficients beta1, beta2 (dimensionless)
        //  beta_{1,2} = (5/96) * (2 -+ sqrt(6/5))
        static const scalar beta1_;
        static const scalar beta2_;

        //- Maximum number of inner Picard iterations on the Larsen
        //  beta cross-coupling per correct() call
        label larsenMaxIters_;

        //- Relative tolerance on the boundary value of phi1_j / phi2_j
        //  for declaring inner Picard iteration converged
        scalar larsenTol_;

        //- Solved auxiliary photon distribution components phi*_{1,j},
        //  one per group, satisfying the first SP3 Helmholtz sub-equation
        PtrList<volScalarField> phi1_j_;

        //- Solved auxiliary photon distribution components phi*_{2,j},
        //  one per group, satisfying the second SP3 Helmholtz sub-equation
        PtrList<volScalarField> phi2_j_;

        //- Partial source fields S_j, one per group, recovered from
        //  phi*_{1,j} and phi*_{2,j} via the SP3 closure combination
        PtrList<volScalarField> Sph_j_;

public:

    //- Runtime type information
    TypeName("threeGroupSP3");

    // Constructors

        //- Construct from mesh and dictionary
        threeGroupSP3(const fvMesh& mesh);

    //- Destructor
    virtual ~threeGroupSP3() = default;


    // Member Functions

        //- Assemble and solve the 6 SP3 sub-equations (2 per group), combine
        //  via the SP3 closure, and sum into the total source field
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

