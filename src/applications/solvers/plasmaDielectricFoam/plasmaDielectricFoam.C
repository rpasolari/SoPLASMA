/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025 Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Application
    plasmaDielectricFoam

Description
    Transient solver for coupled gas (plasma) and dielectric domains developed
    for plasma simulation purposes.

Usage
    \b plasmaDielectricFoam [OPTIONS]

    Example:
        plasmaDielectricFoam -case plasmaDielectricTest

Author
    Rention Pasolari
    Contact: r.pasolari@gmail.com
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "regionProperties.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "fvSolution.H"
#include "solutionControl.H"
#include "mappedPatchBase.H"
#include "pimpleControl.H"

#include "ddSolidSurfaceFluxFvPatchScalarField.H"
#include "foamPlasmaToolkitConstants.H"
#include "plasmaSpecies.H"
#include "plasmaEnergy.H"
#include "plasmaTransport.H"
#include "plasmaTimeControl.H"
#include "plasmaProfiler.H"
// #include "adaptiveFvMesh.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for coupled gas (plasma) and dielectric domains"
        " developed for plasma simulation purposes."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    pimpleControl pimple(gasMesh());
    plasmaTimeControl timeControl(runTime, gasMesh());
    timeControl.setInitialDeltaT(transport);
    
    #include "createCoupledRegions.H"
    // #include "readAMRConfiguration.H"

    #include "readElectricPotentialControls.H"

    #include "reportSimulationSummary.H"
    #include "updateChargeDensity.H"

    // Compute the initial Electric field
    if (poissonSolver == "explicit")
    {
        #include "solveElectricPotential.H"
        #include "calculateElectricField.H"
    }

    runTime.writeNow();

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;
        // #include "syncMultiRegionAMR.H"

        timeControl.adjustDeltaT(transport);

        if (poissonSolver == "semiImplicit")
        {
            while (pimple.loop())
            {
                while (pimple.correct())
                {
                    #include "solveElectricPotential.H"
                    #include "calculateElectricField.H"
                }
                for (nonOrth = 0; nonOrth <= nNonOrthoCorr; ++nonOrth)
                {
                    bool firstIter = (nonOrth == 0);
                    bool finalIter = (nonOrth == nNonOrthoCorr);

                    transport.correct(firstIter, finalIter);
                }

                #include "updateChargeDensity.H"
                #include "updateSurfCharge.H"
            }
        }
        else if (poissonSolver == "explicit")
        {
            plasmaProfiler::start("Solver Sequence", "Explicit Poisson Full Loop");
            int pIter = 0;

            while (pimple.loop())
            {
                pIter++;
                std::string iterLabel = "PIMPLE Iter " + std::to_string(pIter);

                plasmaProfiler::start("Iteration Control", iterLabel);
                gasMesh().update();

                plasmaProfiler::start("Physics", iterLabel + ": Transport");
                for (nonOrth = 0; nonOrth <= nNonOrthoCorr; ++nonOrth)
                {
                    bool firstIter = (nonOrth == 0);
                    bool finalIter = (nonOrth == nNonOrthoCorr);
                    transport.correct(firstIter, finalIter);
                }
                plasmaProfiler::stop("Physics", iterLabel + ": Transport");

                plasmaProfiler::start("Source Terms", iterLabel + ": Volumetric Charge");;
                #include "updateChargeDensity.H"
                plasmaProfiler::stop("Source Terms", iterLabel + ": Volumetric Charge");

                plasmaProfiler::start("Source Terms", iterLabel + ": Surface Charge");
                #include "updateSurfCharge.H"
                plasmaProfiler::stop("Source Terms", iterLabel + ": Surface Charge");

                plasmaProfiler::start("Linear Algebra", iterLabel + ": Poisson Solver");
                #include "solveElectricPotential.H"
                plasmaProfiler::stop("Linear Algebra", iterLabel + ": Poisson Solver");

                plasmaProfiler::start("Physics", iterLabel + ": E-Field Calc");
                #include "calculateElectricField.H"
                plasmaProfiler::stop("Physics", iterLabel + ": E-Field Calc");

                plasmaProfiler::stop("Iteration Control", iterLabel);
            }
            plasmaProfiler::stop("Solver Sequence", "Explicit Poisson Full Loop");
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    plasmaProfiler::report();

    Info<< "End\n" << endl;

    return 0;
}
