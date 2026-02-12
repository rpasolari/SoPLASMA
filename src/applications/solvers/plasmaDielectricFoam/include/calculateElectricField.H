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

E = -fvc::grad(ePotential);
E.correctBoundaryConditions();

forAll(dielectricRegions, i)
{
    volScalarField& ePotentialDielectric = ePotentialDielectricList[i];
    volVectorField& EDielectric = EDielectricList[i];

    EDielectric = -fvc::grad(ePotentialDielectric);
    EDielectric.correctBoundaryConditions();
}
