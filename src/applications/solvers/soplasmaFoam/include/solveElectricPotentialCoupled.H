/*---------------------------------------------------------------------------*\
  File: solveElectricPotentialCoupled.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson/Laplace equations for the electric potential in all
      regions (gas + dielectric), in a coupled way (monolithically). 

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

for (nonOrth = 0; nonOrth <= nonOrthoCorrCoupled; ++nonOrth)
{
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

    if (nonOrth == 0)
    {
        Info<< "Solving for ePotential in coupled regions (monolithic)" << endl;
    }

    fvMatrixAssemblyPtr->solve();

    ePotential.correctBoundaryConditions();

    forAll(dielectricRegions, i)
    {
        ePotentialDielectricList[i].correctBoundaryConditions();
    }

    if (nonOrth == nonOrthoCorrCoupled)
    {
        fvOptions.correct(ePotential);
    }

    fvMatrixAssemblyPtr->clear();
}






