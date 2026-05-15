/*---------------------------------------------------------------------------*\
  File: solveElectricPotentialCoupled.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Solves the Poisson/Laplace equations for the electric potential in all
      regions (fluid + dielectric), in a coupled way (monolithically). 

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

forAll(fluidRegions, i)
{
    ePotentialFluid[i].correctBoundaryConditions();
}
forAll(dielectricRegions, i)
{
    ePotentialDielectric[i].correctBoundaryConditions();
}

for (int nonOrth = 0; nonOrth <= nonOrthoCorrCoupled; ++nonOrth)
{
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

    if (nonOrth == 0)
    {
        Info<< "Solving for ePotential in coupled regions (monolithic)"
            << nl << endl;
    }

    fvMatrixAssemblyPtr->solve();

    fvMatrixAssemblyPtr->clear();
}

Info<< "\n>>> Solved for ePotential in coupled regions" << endl;





