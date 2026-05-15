/*---------------------------------------------------------------------------*\
  File: PoissonEquation.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson equation for the electric potential in the plasma
      (fluid) regions:

          ∇·(ε∇φ) = -ρ_e

      where ε is the permittivity and ρ_e is the charge density. This
      equation couples the plasma charge distribution to the resulting
      electrostatic potential field.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if (!coupled)
{
    Info<< "\nSolving for ePotential in region " 
        << fluidRegions[i].name() << endl;

    for (int nonOrth = 0; nonOrth <= nNonOrthoCorrFluid[i]; ++nonOrth)
    {
        fvScalarMatrix PoissonEquation
        (
            fvm::laplacian(epsilon,ePotential)
            ==
            -eCharge
        );

        auto sp = PoissonEquation.solve();

        if (nonOrth == 0)
        {
            ePotentialNonCoupledResidualFluid[i] = sp.initialResidual();
        }
    }
}
else
{
        fvScalarMatrix PoissonEquation
        (
            fvm::laplacian(epsilon,ePotential)
            ==
            -eCharge
        );

        fvMatrixAssemblyPtr->addFvMatrix(PoissonEquation);
}

