/*---------------------------------------------------------------------------*\
  File: readElectricPotentialControls.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Reads solver and convergence controls for the electric potential
      equation (Poisson/Laplace), for both coupled and segregated modes.

      In segregated (non-coupled) mode, it reads per-region convergence
      tolerances and optional maxIter settings for the outer loop.
      Iterations terminate once all regions fall below tolerance or when maxIter
      is reached.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

List<scalar> ePotentialNonCoupledToleranceFluid(fluidRegions.size());
List<scalar> ePotentialNonCoupledToleranceDielectric(dielectricRegions.size());

List<int> nNonOrthoCorrFluid(fluidRegions.size());
List<int> nNonOrthoCorrDielectric(dielectricRegions.size());

int nonOrthoCorrCoupled = 0;
scalar ePotentialNonCoupledLoopMaxIter = 0;
bool foundAnyRegionMaxIter = false;

forAll(fluidRegions, i)
{
    const fvMesh& mesh = fluidRegions[i];

    const dictionary epCtrl =
        mesh.solutionDict().subDict("ePotentialControls");

    // Always read — used by both coupled and segregated solve
    nNonOrthoCorrFluid[i] =
        epCtrl.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    nonOrthoCorrCoupled = max(nonOrthoCorrCoupled, nNonOrthoCorrFluid[i]);

    // Only used in segregated mode
    const dictionary nonCoupledCtrl =
        epCtrl.subOrEmptyDict("nonCoupledResidualControl");

    ePotentialNonCoupledToleranceFluid[i] =
        nonCoupledCtrl.getOrDefault<scalar>("tolerance", 1.0e-6);

    if (nonCoupledCtrl.found("maxIter"))
    {
        foundAnyRegionMaxIter = true;
        scalar regionMax =
            nonCoupledCtrl.get<scalar>("maxIter");
        ePotentialNonCoupledLoopMaxIter =
            max(ePotentialNonCoupledLoopMaxIter, regionMax);
    }
}

forAll(dielectricRegions, i)
{
    const fvMesh& mesh = dielectricRegions[i];

    const dictionary epCtrl =
        mesh.solutionDict().subDict("ePotentialControls");

    nNonOrthoCorrDielectric[i] =
        epCtrl.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    nonOrthoCorrCoupled = max(nonOrthoCorrCoupled, nNonOrthoCorrDielectric[i]);

    const dictionary nonCoupledCtrl =
        epCtrl.subOrEmptyDict("nonCoupledResidualControl");

    ePotentialNonCoupledToleranceDielectric[i] =
        nonCoupledCtrl.getOrDefault<scalar>("tolerance", 1.0e-6);

    if (nonCoupledCtrl.found("maxIter"))
    {
        foundAnyRegionMaxIter = true;
        scalar regionMax =
            nonCoupledCtrl.get<scalar>("maxIter");
        ePotentialNonCoupledLoopMaxIter =
            max(ePotentialNonCoupledLoopMaxIter, regionMax);
    }
}

// Default outer iterations if none specified
if (!foundAnyRegionMaxIter)
{
    ePotentialNonCoupledLoopMaxIter = 100;
}

// Print settings only in segregated mode
if (!coupled)
{
    Info<< "*** Segregated electric potential controls:" << nl
        << "   - Max iterations: " << ePotentialNonCoupledLoopMaxIter << nl;

    forAll(fluidRegions, i)
    {
        Info<< "   - (Fluid region) " << fluidRegions[i].name()
            << " : nNonOrth=" << nNonOrthoCorrFluid[i]
            << " , tolerance: " << ePotentialNonCoupledToleranceFluid[i] << nl;
    }
    forAll(dielectricRegions, i)
    {
        Info<< "   - (Dielectric region) " << dielectricRegions[i].name()
            << " : nNonOrth=" << nNonOrthoCorrDielectric[i]
            << " , tolerance: "
            << ePotentialNonCoupledToleranceDielectric[i] << nl;
    }
    Info<< endl;
}
else
{
    Info<< "*** Monolithic electric potential controls:" << nl
        << "   - Using maximum nNonOrthogonalCorrectors across all regions: "
        << nonOrthoCorrCoupled << nl
        << endl;
}
