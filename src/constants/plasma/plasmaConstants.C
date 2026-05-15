/*---------------------------------------------------------------------------*\
  File: plasmaConstants.C
  Part of: SoPLASMA

  Description:
      Implementation file for plasma-related physical constants declared in
      plasmaConstants.H.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaConstants.H"
#include "mathematicalConstants.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::cLight,
    dimensionedScalar
    (
        "cLight",
        dimensionSet(0, 1, -1, 0, 0, 0, 0),
        2.99792458e8
    ),
    constantplasmacLight,
    "cLight"
);

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::eCharge,
    dimensionedScalar
    (
        "eCharge",
        dimensionSet(0, 0, 1, 0, 0, 1, 0),
        1.602176634e-19
    ),
    constantplasmaeCharge,
    "eCharge"
);

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::eMass,
    dimensionedScalar
    (
        "eMass",
        dimensionSet(1, 0, 0, 0, 0, 0, 0),
        9.1093837139e-31
    ),
    constantplasmaeMass,
    "eMass"
);

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::mu0,
    dimensionedScalar
    (
        "mu0",
        dimensionSet(1, 1, -2, 0, 0, -2, 0),
        4.0*mathematical::pi*1e-07
    ),
    constantplasmamu0,
    "mu0"
);

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::epsilon0,

    dimensionedScalar
    (
        "epsilon0",
        dimensionedScalar
        (
            "epsilon0",
            dimensionSet(0, 0, 0, 0, 0),
            1.0
        )
       /(plasma::mu0*sqr(plasma::cLight))
    ),
    constantplasmaepsilon0,
    "epsilon0"
);

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::kappaCoulomb,

    dimensionedScalar
    (
        "kappaCoulomb",
        dimensionedScalar
        (
            "kappaCoulomb",
            dimensionSet(0, 0, 0, 0, 0),
            1.0/(4.0*mathematical::pi)
        )
       /plasma::epsilon0
    ),

    constantplasmakappaCoulomb,
    "kappaCoulomb"
);

defineDimensionedConstantWithDefault
(
    plasma::group,
    plasma::kappaBoltzmann,

    dimensionedScalar
    (
        "kappaBoltzmann",
        dimensionedScalar
        (
            "kappaBoltzmann",
            dimensionSet(1, 2, -2, -1, 0, 0, 0),
            1.380649e-23
        )
    ),

    constantplasmakappaBoltzmann,
    "kappaBoltzmann"
);

} // End namespace constant
} // End namespace Foam

// ************************************************************************* //
