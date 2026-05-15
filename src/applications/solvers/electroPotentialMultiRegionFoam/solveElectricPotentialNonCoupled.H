/*---------------------------------------------------------------------------*\
  File: solveElectricPotentialNonCoupled.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson/Laplace equations for the electric potential in all
      regions (fluid + dielectric), in a non-coupled way (segregate).

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info<< "Solving for ePotential in " 
    << "non-coupled regions (segregate)" << endl;

for(int outer = 1; outer <= ePotentialNonCoupledLoopMaxIter; ++outer)
{

    Info<< "\n  -ePotential iteration (" << outer << "/" 
        << ePotentialNonCoupledLoopMaxIter << ")" << endl;

    // Reset residual lists each iteration
    forAll(ePotentialNonCoupledResidualFluid, i)
        ePotentialNonCoupledResidualFluid[i] = Foam::GREAT;

    forAll(ePotentialNonCoupledResidualDielectric, i)
        ePotentialNonCoupledResidualDielectric[i] = Foam::GREAT;

    forAll(fluidRegions, i)
    {
        fvMesh& mesh = fluidRegions[i];
        #include "setRegionFluidFields.H"
        #include "PoissonEquation.H"
    }

    forAll(dielectricRegions, i)
    {
        fvMesh& mesh = dielectricRegions[i];
        #include "setRegionDielectricFields.H"
        #include "LaplaceEquation.H"
    }

    // Check convergence
    bool allOK = true;

    forAll(ePotentialNonCoupledResidualFluid, i)
    {
        if
        (
            ePotentialNonCoupledResidualFluid[i]
            >= ePotentialNonCoupledToleranceFluid[i]
        )
        {
            allOK = false;
        }
    }

    forAll(ePotentialNonCoupledResidualDielectric, i)
    {
        if
        (
            ePotentialNonCoupledResidualDielectric[i]
            >= ePotentialNonCoupledToleranceDielectric[i]
        )
        {
            allOK = false;
        }
    }

    if (allOK)
    {
        Info<< "\n>>> ePotential converged in "
            << outer << " iterations." << nl;
        break;
    }
    else if (outer == ePotentialNonCoupledLoopMaxIter)
    {
        Info<< "\n>>> WARNING: ePotential did NOT converge after "
            << outer << " iterations." << nl;
    }
}



