/*---------------------------------------------------------------------------*\
  File: solveElectricPotentialNonCoupled.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson/Laplace equations for the electric potential in all
      regions (gas + dielectric), in a non-coupled way (segregate).

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

if (nTotalRegions > 1)
{
    Info<< "Solving for ePotential in " 
        << "non-coupled regions (segregate)" << endl;
}

for(int outer = 1; outer <= ePotentialNonCoupledLoopMaxIter; ++outer)
{

    if (ePotentialNonCoupledLoopMaxIter > 1)
    {
        Info<< "  -ePotential iteration (" << outer << "/" 
            << ePotentialNonCoupledLoopMaxIter << ")" << endl;
    }

    // Reset residual lists each iteration
    ePotentialNonCoupledResidual = Foam::GREAT;

    forAll(ePotentialNonCoupledResidualDielectricList, i)
        ePotentialNonCoupledResidualDielectricList[i] = Foam::GREAT;

    {
        fvMesh& mesh = gasMesh();
        #include "PoissonEquation.H"
    }

    forAll(dielectricRegions, i)
    {
        fvMesh& mesh = dielectricRegions[i];
        #include "setRegionDielectricFields.H"
        #include "LaplaceEquation.H"
    }

    if (nTotalRegions == 1) break;

    // Check convergence
    bool allOK = (ePotentialNonCoupledResidual < ePotentialNonCoupledTolerance);

    forAll(ePotentialNonCoupledResidualDielectricList, i)
    {
        if
        (
            ePotentialNonCoupledResidualDielectricList[i]
            >= ePotentialNonCoupledToleranceDielectricList[i]
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



