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

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if (!coupled)
{
    Info<< "Solving for ePotential in region " 
        << dielectricRegions[i].name() << endl;

    for (nonOrth = 0; nonOrth <= nNonOrthoCorrDielectricList[i]; ++nonOrth)
    {
        fvScalarMatrix LaplaceEquation
        (
            fvm::laplacian(epsilonDielectric,ePotentialDielectric)
        );

        if (nonOrth < nonOrthoCorrCoupled)
        {
            LaplaceEquation.relax();
        }

        auto sp = LaplaceEquation.solve();

        if (nonOrth == 0)
        {
            ePotentialNonCoupledResidualDielectricList[i] 
                = sp.initialResidual();
        }
    }
}
else
{
        fvScalarMatrix LaplaceEquation
        (
            fvm::laplacian(epsilonDielectric,ePotentialDielectric)
        );

        if (nonOrth < nonOrthoCorrCoupled)
        {
            LaplaceEquation.relax();
        }
        
        fvMatrixAssemblyPtr->addFvMatrix(LaplaceEquation);
}
