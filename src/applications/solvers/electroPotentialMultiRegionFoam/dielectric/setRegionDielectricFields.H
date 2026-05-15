/*---------------------------------------------------------------------------*\
  File: setRegionSolidFields.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from setRegionSolidFields.H in OpenFOAM v2412.
      Defines convenient field references for each dielectric region.
      Associates region-specific variables with their corresponding
      data structures created in createDielectricFields.H.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

volScalarField& ePotential = ePotentialDielectric[i];
volScalarField* surfCharge = surfChargeDielectric.set(i)
                                ? &surfChargeDielectric[i]
                                : nullptr;

dimensionedScalar& epsilon = epsilonDielectric[i];



