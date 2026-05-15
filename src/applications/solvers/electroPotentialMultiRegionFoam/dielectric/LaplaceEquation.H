/*---------------------------------------------------------------------------*\
  File: LaplaceEquation.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Laplace equation for the electric potential in the dielectric
      regions:

          ∇·(ε∇φ) = 0

      where ε is the permittivity. This represents electrostatic equilibrium
      in a charge-free region (ρ_e = 0). 

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if (!coupled)
{
    Info<< "\nSolving for ePotential in region " 
        << dielectricRegions[i].name() << endl;

    for (int nonOrth = 0; nonOrth <= nNonOrthoCorrDielectric[i]; ++nonOrth)
    {
        fvScalarMatrix LaplaceEquation
        (
            fvm::laplacian(epsilon,ePotential)
        );

        auto sp = LaplaceEquation.solve();

        if (nonOrth == 0)
        {
            ePotentialNonCoupledResidualDielectric[i] 
                = sp.initialResidual();
        }
    }
}
else
{
        fvScalarMatrix LaplaceEquation
        (
            fvm::laplacian(epsilon,ePotential)
        );

        fvMatrixAssemblyPtr->addFvMatrix(LaplaceEquation);
}
