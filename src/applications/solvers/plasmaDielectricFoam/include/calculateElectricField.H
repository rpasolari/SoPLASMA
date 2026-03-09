/*---------------------------------------------------------------------------*\
  File: calculateElectricField.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Description:
    Calculates the electric field from the electric potential, for both the gas
    and the dielectric regions.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

plasmaProfiler::start("Physics : E-Field Grad (Gas)");
E = -fvc::grad(ePotential);
plasmaProfiler::stop("Physics : E-Field Grad (Gas)");

// --- 2. Gas Phase Boundary Correction ---
plasmaProfiler::start("Physics : E-Field BC Correct (Gas)");
E.correctBoundaryConditions();
plasmaProfiler::stop("Physics : E-Field BC Correct (Gas)");

plasmaProfiler::start("Physics : All Dielectrics (Grad + BC)");

forAll(dielectricRegions, i)
{
    volScalarField& ePotentialDielectric = ePotentialDielectricList[i];
    volVectorField& EDielectric = EDielectricList[i];

    EDielectric = -fvc::grad(ePotentialDielectric);
    EDielectric.correctBoundaryConditions();
}

plasmaProfiler::stop("Physics : All Dielectrics (Grad + BC)");
