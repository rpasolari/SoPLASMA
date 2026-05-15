/*---------------------------------------------------------------------------*\
  File: plasmaConstants.H
  Part of: SoPLASMA

  Description:
      Header providing plasma-related physical constants for use throughout
      the SoPLASMA. These constants follow the OpenFOAM conventions
      for dimensioned types and are intended for shared use by solvers,
      models, and boundary conditions.

  Namespace:
      Foam::constant::plasma

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef plasmaConstants_H
#define plasmaConstants_H

#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{
namespace plasma
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Group name for plasma constants
    constexpr const char* const group = "plasma";

    //- Speed of light in vacuum: default SI units: [m/s]
    extern const dimensionedScalar cLight;    

    //- Elementary charge: default SI units: [C]
    extern const dimensionedScalar eCharge;

    //- Electron mass: default SI units: [kg]
    extern const dimensionedScalar eMass;

    //- Magnetic constant/permeability of free space: default SI units: [H/m]
    extern const dimensionedScalar mu0;

    //- Dielectric permittivity of free space: default SI units: [F/m]
    extern const dimensionedScalar epsilon0;

    //- Coulomb constant: default SI units: [N.m2/C2]
    extern const dimensionedScalar kappaCoulomb;

    //- Coulomb constant: default SI units: [N.m2/C2]
    extern const dimensionedScalar kappaBoltzmann;

    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace plasma
} // End namespace constant
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


