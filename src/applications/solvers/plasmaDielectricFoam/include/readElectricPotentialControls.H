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

scalar ePotentialNonCoupledTolerance;
int nNonOrthoCorr;
int nonOrth;

List<scalar> ePotentialNonCoupledToleranceDielectricList
(
    dielectricRegions.size()
);
List<int> nNonOrthoCorrDielectricList(dielectricRegions.size());

int nonOrthoCorrCoupled = 0;
scalar ePotentialNonCoupledLoopMaxIter = 0;
bool foundAnyRegionMaxIter = false;

{
    const fvMesh& mesh = gasMesh();

    const dictionary epCtrl =
        mesh.solutionDict().subDict("ePotentialControls");

    // Always read - used by both coupled and segregated solve
    nNonOrthoCorr =
        epCtrl.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    nonOrthoCorrCoupled = max(nonOrthoCorrCoupled, nNonOrthoCorr);

    // Only used in segregated mode
    const dictionary nonCoupledCtrl =
        epCtrl.subOrEmptyDict("nonCoupledResidualControl");

    ePotentialNonCoupledTolerance =
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

    nNonOrthoCorrDielectricList[i] =
        epCtrl.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    nonOrthoCorrCoupled = max
    (
        nonOrthoCorrCoupled,
        nNonOrthoCorrDielectricList[i]
    );

    const dictionary nonCoupledCtrl =
        epCtrl.subOrEmptyDict("nonCoupledResidualControl");

    ePotentialNonCoupledToleranceDielectricList[i] =
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

if (nTotalRegions == 1)
{
    ePotentialNonCoupledLoopMaxIter = 1;
}
