/*---------------------------------------------------------------------------*\
  File: PoissonEquation.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson equation for the electric potential in the plasma
      (gas) region:

          ∇·(ε∇φ) = -ρ_e

      where ε is the permittivity and ρ_e is the charge density. This
      equation couples the plasma charge distribution to the resulting
      electrostatic potential field.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

const dimensionedScalar deltaT = gasMesh().time().deltaT();

if (!coupled)
{
    Info<< "Solving for ePotential in region " 
        << gasRegionName << endl;

    for (nonOrth = 0; nonOrth <= nNonOrthoCorr; ++nonOrth)
    {
        if (poissonSolver == "explicit")
        {
            fvScalarMatrix PoissonEquation
            (
                fvm::laplacian(epsilon,ePotential)
                ==
                -chargeDensity
            );

            fvOptions.constrain(PoissonEquation);

            if (nonOrth < nNonOrthoCorr)
            {
                PoissonEquation.relax();
            }

            auto sp = PoissonEquation.solve();

            if (nonOrth == 0)
            {
                ePotentialNonCoupledResidual = sp.initialResidual();
            }
        }
        else if (poissonSolver == "semiImplicit")
        {
            fvScalarMatrix PoissonEquation(
                fvm::laplacian(
                    epsilon + deltaT * transport.electricalConductivity(),
                    ePotential
                )
            ==
                -chargeDensity
                - deltaT * transport.diffusiveChargeSource()
            );
            
            fvOptions.constrain(PoissonEquation);

            if (nonOrth < nNonOrthoCorr)
            {
                PoissonEquation.relax();
            }

            auto sp = PoissonEquation.solve();

            if (nonOrth == 0)
            {
               ePotentialNonCoupledResidual = sp.initialResidual();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unknown Poisson solver: " << poissonSolver
                << exit(FatalError);
        }
    }

    fvOptions.correct(ePotential);
}
else
{
    if (poissonSolver == "explicit")
    {
        fvScalarMatrix PoissonEquation
        (
            fvm::laplacian(epsilon,ePotential)
            ==
            -chargeDensity
        );

        fvOptions.constrain(PoissonEquation);

        if (nonOrth < nonOrthoCorrCoupled)
        {
            PoissonEquation.relax();
        }

        fvMatrixAssemblyPtr->addFvMatrix(PoissonEquation);
    }
    else if (poissonSolver == "semiImplicit")
    {
        fvScalarMatrix PoissonEquation(
            fvm::laplacian(
                epsilon + deltaT * transport.electricalConductivity(),
                ePotential
            )
        ==
            -chargeDensity
            -deltaT * transport.diffusiveChargeSource()
        );

        fvOptions.constrain(PoissonEquation);

        if (nonOrth < nonOrthoCorrCoupled)
        {
            PoissonEquation.relax();
        }

        fvMatrixAssemblyPtr->addFvMatrix(PoissonEquation);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown Poisson solver: " << poissonSolver
            << exit(FatalError);
    }
}
